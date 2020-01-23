# -*- coding: utf-8 -*-
=begin
= 無次元化した Boussinesq, 鉛直ー水平２次元モデル
x: cyclic, 範囲 0..2pi
z: bounded, 範囲 0..pi

== 支配方程式
  ζ_t + α(u ζ_x+w ζ_z) = -αb_x - ν((-Δ)^nl)ζ
  b_t + α(u b_x+w b_z) + αN2 w = αQ - κ((-Δ)^nl)ζ
  ζ = Δφ = α^2 φ_xx + φ_zz,  u=φ_z, w=-φ_x


=end

require "narray"
require "numru/fftw3"
require "numru/gphys"
require File.dirname(File.expand_path(__FILE__)) + "/fftw_ext.rb"

module NumRu
  class Boussinesq2D
    attr_reader :b, :Q, :dx, :dz, :t, :dimensional_out, :n2, :n2e

    # Boussineq では nondim_H はダミー(Anelasticと引数を合わせて切り替えやすく)
    def initialize(nx, nz, a, n2, nu, kp, nl, linear: false, hydro: false,
                   hs: nil, ts: nil, w_hips: false, k_v_h: 1.0, nondim_H: nil)
      @nx = nx
      @nz = nz        # 鉛直範囲分割数
      @nz1 = nz - 1   # モデル中のデータの鉛直格子点数＝両端を除く
      @nze = nz + 1   # 両端を加えた場合の鉛直格子点数(GPhys/ファイル出力用)
      @a = a     # alpha = H/L
      if n2.respond_to?(:length) && n2.length != @nze
        raise("invalid n2 length #{n2.length} != #{@nze}")
      end
      @n2e = n2   # 無次元N2 (定数または長さ nz+1 のNArray: 第0..nz層まで)
      if n2.respond_to?(:length)
        @n2 = @n2e[false,1...@nz]
      else
        @n2 = @n2e
      end      
      @nu = nu   # 無次元動粘性係数 (can be hyper, 最大波数の減衰率)
      @kp = kp   # 無次元熱伝導係数 (can be hyper, 最大波数の減衰率)
      @nl = nl   # 粘性・伝導のラプラシアンの次数 (1が普通だが，>1 で hyper)
      @linear = linear # true なら移流項無視
      @hydro = hydro  # trueなら静水圧(laplacianでm=0; linear用)
      @hs = hs  # 鉛直スケールパラメタ Z (座標軸設定用) [m]  (hs, tsは同時指定)
      @hs_km = @hs/1000.0 if @hs # [km] バージョン
      @ts = ts  # 時間スケールパラメタ 1/N0 (座標軸設定用) [s](hs, tsは同時指定)
      if @hs && @ts.nil? || @ts && @hs.nil?
        raise("hs and ts must be specified together")
      end
      @dimensional_out = !@hs.nil?
      
      phi = NArray.float(nx,@nz1)   # stream function (本来内部変数でなくていい)
      b = NArray.float(nx,@nz1)     # buoyancy (本来内部変数でなくていい)
      @Q = NArray.float(nx,@nz1)     # heating (set_Qで設定)
      @step = 0
      @t = 0.0
      @phi_s = g2s(phi)   # stream function
      @b_s = g2s(b)       # buoyancy
      @w_hips = w_hips
      @k_v_h = k_v_h  # 鉛直拡散粘性を調整 kp_v = k_v_h * kp, nu_v = k_v_h * nu
      
      @k = NArray.float(nx).indgen!
      @k[nx/2+1..-1] = nx - @k[nx/2+1..-1]
      @m = NArray.float(@nz1).indgen! + 1
      @m_d = @m.newdim(0)
      @k_dc = -NArray.float(nx/2-1).indgen!-1   # -1,-2,...,-(nx/2-1)
      @k_ds = -@k_dc[-1..0]    # -(nx/2-1),...-2,-1
      @k_dc = @k_dc.newdim!(1)
      @k_ds = @k_ds.newdim!(1)
      #p @k_dc, @k_ds
      if @hydro
        @lapla = -( @k.newdim(-1)*0 + @m_d**2 )
      else
        @lapla = -( @a**2 * @k.newdim(-1)**2 + @m_d**2 )
      end
      lap_dif = -( @a**2 * @k.newdim(-1)**2 + @k_v_h * @m.newdim(0)**2 )**nl
      @zeta_diff = (@nu/lap_dif.min.abs) * lap_dif
      @b_diff = (@kp/lap_dif.min.abs) * lap_dif
      @ilapla = 1.0/@lapla
      # @ilapla[ @lapla.eq(0.0) ] = 0  # 存在しないはず
      # p "$$$", @nz, lap_dif  #, @ilapla, @ilapla.max, @ilapla.min
      
      @dx = 2*Math::PI/nx   # 無次元
      @dz = Math::PI/@nz     # 無次元

      if @w_hips  # heat-induced passive scalar
        @hips_s = @b_s.clone
      end

      @x = ( NArray.float(nx).indgen! - nx/2 )* @dx
      @z = ( NArray.float(@nz1).indgen! + 1 ) * @dz
      @ze = NArray.float(@nze).indgen! * @dz

        # < GPhys 化関連 >

      if @n2e.is_a?(Numeric)
        @n2_int = (@n2e*@ze).newdim!(0)    # @n2_int : n2の積分量 (「温位」用)
      else
        n2_int = NArray.float(@nze)
        n2mid = (@n2e[1..-1] + @n2e[0..-2])
        for k in 1...nz
          n2_int[k] = n2_int[k-1] + n2mid[k-1] * @dz   # 第2層以上は台形公式で設定
        end
        @n2_int = n2_int.newdim!(0)
      end
      @Qupdater = nil   # Proc にすれば加熱 @Q は状態依存＝時間積分時に更新

      if @dimensional_out  # 次元付きで
        @x_d = @x * (@hs_km/@a)
        @ze_d = @ze * @hs_km
        vx = VArray.new(@x_d,{"long_name"=>"x","units"=>"km",
                             "modulo"=>[@nx*@dx*(@hs_km/@a)]},"x")
        vz = VArray.new(@ze_d,{"long_name"=>"z","units"=>"km"},"z")
      else  # 無次元で
        vx = VArray.new(@x,{"long_name"=>"x","units"=>""},"x")
        vz = VArray.new(@ze,{"long_name"=>"z","units"=>""},"z")
      end

      ax = Axis.new.set_pos(vx)
      az = Axis.new.set_pos(vz)
      @grid = Grid.new(ax,az)
      @gphi = GPhys.new(@grid,
                        VArray.new(g_ze(phi),{"long_name"=>"phi","units"=>""},"phi"))
      @gb = GPhys.new(@grid, VArray.new(g_ze(b),{"long_name"=>"b","units"=>""},"b") )

      if @dimensional_out
        @gphi.units = "m2.s-1"
        @gb.units = "m.s-2"
      end

      if @w_hips  # heat-induced passive scalar (常に無次元)
        @ghips = GPhys.new(@grid,
             VArray.new(g_ze(b),{"long_name"=>"passive scalar","units"=>""},"hips") )
      end

    end

    # 時間 (@ts が定義されていれば次元付き時間）
    def t_d
      if @ts
        @t * @ts
      else
        @t
      end
    end

    # t_d の単位
    def t_d_un
      if @ts
        "s"
      else
        ""
      end
    end

    def phi
      s2g(@phi_s)
    end
    def b
      s2g(@b_s)
    end
    def gphi
      phi = s2g(@phi_s)
      phi *= (@hs**2/@ts) if @dimensional_out
      @gphi[true,1..@nz1] = phi
      @gphi
    end
    def gb
      b = s2g(@b_s)
      b *= (@hs/@ts**2) if @dimensional_out
      @gb[true,1..@nz1] = b
      @gb
    end
    def gbb  # 温位相当量 (b + \int_z αN2) : 浮力ベース，Lagrange保存
      if @dimensional_out
        gbb = gb + @n2_int * (@hs/@ts**2)
      else
        gbb = gb + @n2_int
      end
      gbb.name = "theta"
      gbb.long_name = "Boussinesq pot temp"
      gbb
    end
    def ghips
      raise("Not initialized with heat-induced passive scalar is") if !@w_hips
      hips = s2ge(@hips_s)
      @ghips.replace_val(hips)
      @ghips
    end

    def gzeta
      zeta = s2ge(@phi_s*@lapla)
      if @dimensional_out
        zeta /= @ts
        un = "s-1"
      else
        un = ""
      end
      gzeta = GPhys.new(@grid,
                    VArray.new(zeta,{"long_name"=>"zeta","units"=>un},"zeta"))
      gzeta
    end

    def gu
      u = phi_s2u_e(@phi_s)
      if @dimensional_out
        u *= (@hs/@ts)
        un = "m/s"
      else
        un = ""
      end
      gu = GPhys.new(@grid,
                    VArray.new(u,{"long_name"=>"U","units"=>un},"u"))
      gu
    end

    def gw
      w = g_ze( phi_s2w(@phi_s) )
      if @dimensional_out
        w *= (@a*@hs/@ts)
        un = "m/s"
      else
        un = ""
      end
      gw = GPhys.new(@grid,
                    VArray.new(w,{"long_name"=>"W","units"=>un},"w"))
      gw
    end

    def gq
      if @dimensional_out
        q = @Q * (@a*@hs/@ts**3)
        un = "m.s-3"
      else
        q = @Q.clone
        un = ""
      end
      q = g_ze( q )
      gq = GPhys.new(@grid,
              VArray.new(q,{"long_name"=>"diabatic heating","units"=>un},"Q"))
      gq
    end

    # grid to spectral data
    def g2s(g)
      s = FFTW3.sin_fw(g,1)
      s = FFTW3.fft_fw(s,0)
      s
    end

    # spectral to grid data
    def s2g(s)
      g = FFTW3.fft_bk(s,0)
      g = FFTW3.sin_bk(g,1)
      g
    end

    def g_ze(g)
      ge = NArray.float(@nx,@nze)
      ge[true,1...@nz] = g
      ge
    end

    def g_zc(ge)
      ge[true,1...@nz]
    end

    # spectral to grid data (extended)
    def s2ge(s)
      g_ze( s2g(s) )
    end

    def set_sponge(nl, rate)
      @sponge = true
      fact = NArray.float(1,@nz1)
      fact[0,-nl..-1] = (NArray.float(nl).indgen!+1) * (-rate/nl.to_f)
      @tend_sponge = Proc.new{|f|
        f * fact
      }
    end
    
    def set_Q(q=nil, amp=1.0)
      if q.nil?
        nw = 5
        @Q = NArray.float(@nx,@nz1)
        xr = ( (NArray.float(2*nw+1).indgen-nw)*(Math::PI/2/nw) ).newdim(-1)
        @Q[ @nx/2-1-nw..@nx/2-1+nw, true ] = NMath.cos( xr ) * amp
        @Q = @Q * NMath.sin( @z.newdim(0) )
        #@Q[true,0] = 0   # 丸めで完全にはゼロにならいのを強制的にゼロに
        #@Q[true,-1] = 0
      else
        @Q = q
      end
    end

    def clear_Q
      @Q.fill!(0)
    end

    # amp2 は鉛直波数２の振幅：(ampが正のときに) 正->bottom heavy, 負->top heavy
    def set_Q_2(n_ul,nw,amp, amp2=nil)
      q = NArray.float(@nx,@nz1)
      nw1 = nw-1
      xr = ( (NArray.float(2*nw-1).indgen-nw1)*(Math::PI/2/nw) ).newdim(-1)
      ns = [n_ul,@nz1].min
      zr = (NArray.float(ns).indgen!+1).newdim(0) * (Math::PI/n_ul)
      q[@nx/2-nw1..@nx/2+nw1, 0...ns] = amp * NMath.cos(xr) * NMath.sin(zr)
      if amp2
        q[@nx/2-nw1..@nx/2+nw1, 0...ns] += amp2*NMath.cos(xr) * NMath.sin(2*zr)
      end
      @Q = q
    end

    def set_Q_3(n_ul,nw,amp)
      q = NArray.float(@nx,@nz1)
      nw1 = nw-1
      ns = [n_ul,@nz1].min
      xr = ( (NArray.float(2*nw-1).indgen-nw1)*(Math::PI/2/nw) ).newdim(-1)
      zr = (NArray.float(ns).indgen!+1).newdim(0) * (Math::PI / n_ul)
      zshape = NMath.sin(zr) * (1-zr/Math::PI)  # z とともに線形に減らす
                                                # 上層での潜熱放出難に対応
      q[@nx/2-nw1..@nx/2+nw1, 0...ns] = amp * NMath.cos(xr) * zshape
      @Q = q
    end

    # 簡易パラメタリゼーション版
    # 中心付近の領域に N^2 w に比例する加熱。最低値を定め静止状態から誘起する。
    #
    # * boostQamp  最低加熱のピーク値 (正値を指定。boost_untilを指定すれば，その時間まで)
    #   (ただし，w>0の場合のみ加熱：冷却はしない)
    def set_Q_4_p(n_ul, nw, fact: 1.1, boostQamp: 0.1, boost_until: nil)
      @Q = NArray.float(@nx,@nz1)
      nw1 = nw-1
      xr = ( (NArray.float(2*nw-1).indgen-nw1)*(Math::PI/2/nw) ).newdim(-1)
      ns = [n_ul,@nz1].min
      zr = (NArray.float(ns).indgen!+1).newdim(0) * (Math::PI / n_ul)
      @qxran = @nx/2-nw1..@nx/2+nw1
      @qzran = 0...ns
      @Qmin = boostQamp * NMath.cos(xr) * NMath.sin(zr)
      ##@n_lt_thres = 0   # メモ的にカウントするため
      @Qupdater_with_w = Proc.new{|w|
        #w = -s2g_del_x(@phi_s)
        q = fact * (@n2*w)[@qxran, @qzran]
        if boost_until.nil? || @t < boost_until  # 値が @Qmin を下回らないようにする
          mask = q.lt(@Qmin)
          n_lt_thres = mask.count_true
          if n_lt_thres > 0
            q[mask] = @Qmin[mask]
          end
        else
          mask = q.lt(0)
          q[mask] = 0 if mask.count_true > 0
        end
        @Q[@qxran, @qzran] = q
      }
    end

    def update_Q(w: nil)
      if @Qupdater_with_w
        w = -s2g_del_x(@phi_s) unless w
        @Qupdater_with_w.call(w)
        #printf "## update_Q  %f  %f  %d\n", @Q[@qxran, @qzran].max,
        #       w.max, @n_lt_thres
      end
    end
    
    # スペクトルデータを x で微分し格子点に
    def s2g_del_x(s)
      s = s_del_x(s)
      g = FFTW3.fft_bk(s,0)
      g = FFTW3.sin_bk(g,1)
      g
    end

    # スペクトルデータを z で微分し格子点に (上下に伸びるので必要に応じg_zcを適用せよ)
    def s2g_del_z(s)
      s = s * @m_d
      g = FFTW3.fft_bk(s,0)
      g = FFTW3.cos_bk_sinout(g,1)
      g
    end

    # スペクトルのまま x で微分
    def s_del_x(s)
      sc = s[@nx/2-1..1,true]       #  cos部をひっくり返して
      ss = s[@nx-1..@nx/2+1,true]   # -sin部をひっくり返して
      sk = NArray.float(@nx,@nz1)
      sk[1..@nx/2-1,true] = ss*@k_dc       #  cos部 := -sin部 * (-k)
      sk[@nx/2+1..@nx-1,true] = sc*@k_ds   # -sin部 := cos部 * k
      sk
    end

    # スペクトルのままラプラシアン
    def s_lapla(s)
      s = s * @lapla
    end

    def phi_s2u_e(phi_s)
      u = s2g_del_z(phi_s)            # Anelastic の場合に異なる
    end
    def phi_s2u(phi_s)
      u = g_zc( s2g_del_z(phi_s) )    # Anelastic の場合に異なる
    end
    def phi_s2w(phi_s)
      w = -s2g_del_x(phi_s)
    end
    def zeta_z_s_4adv(zeta_s) # ζ_t の移流項における w にかかるもの
      g_zc(s2g_del_z(zeta_s))    # ここではζ_z. Anelastic の場合に異なる
    end
    def adv_factor
      @a
    end

    # tendency
    def tend(phi_s=nil, b_s=nil, hips_s=nil)
      phi_s = @phi_s if phi_s.nil?
      b_s = @b_s if b_s.nil?
      u = phi_s2u(phi_s)
      w = phi_s2w(phi_s)
      zeta_s = s_lapla(phi_s)
      update_Q(w: w)
      af = adv_factor
      if @linear
        #< linear (移流項なし) >
        zeta_t = (-@a)*s_del_x(b_s) + zeta_s * @zeta_diff
        b_t = g2s( - @a*@n2*w + @a*@Q ) + b_s * @b_diff
      else
        #< nonlinear (w/ full advection) >
        adv_zeta = (-u * s2g_del_x(zeta_s) -w * zeta_z_s_4adv(zeta_s) ) * af
        adv_b = (-u * s2g_del_x(b_s) -w * g_zc(s2g_del_z(b_s))) * af
        zeta_t = g2s(adv_zeta) + (-@a)*s_del_x(b_s) + zeta_s * @zeta_diff
        b_t = g2s( adv_b - @a*@n2*w + @a*@Q ) + b_s * @b_diff
      end
      phi_t = zeta_t * @ilapla
      if @sponge
        phi_t += g2s( @tend_sponge.call(s2g(phi_s)) )
        b_t += g2s( @tend_sponge.call(s2g(b_s)) )
      end
      tend = [ phi_t, b_t ]
      if @w_hips
        hips_s = @hips_s if hips_s.nil?
        adv_h = (-u * s2g_del_x(hips_s) -w * g_zc(s2g_del_z(hips_s))) * af
        hips_t = g2s( adv_h + @a*@Q ) + hips_s * @b_diff
                # bとの違いはN2 wがないことと，linearでも移流すること
        tend.push( hips_t )
      end
      tend
    end

    def integ_RK4(dt)
      dt2 = dt/2
      phi_t1, b_t1, hips_t1 = tend
      if @w_hips
        phi_t2, b_t2, hips_t2 = tend(@phi_s+phi_t1*dt2, @b_s+b_t1*dt2,
                                     @hips_s+hips_t1*dt2)
        phi_t3, b_t3, hips_t3 = tend(@phi_s+phi_t2*dt2, @b_s+b_t2*dt2,
                                     @hips_s+hips_t2*dt2)
        phi_t4, b_t4, hips_t4 = tend(@phi_s+phi_t3*dt, @b_s+b_t3*dt,
                                     @hips_s+hips_t3*dt)
      else
        phi_t2, b_t2 = tend(@phi_s+phi_t1*dt2, @b_s+b_t1*dt2)
        phi_t3, b_t3 = tend(@phi_s+phi_t2*dt2, @b_s+b_t2*dt2)
        phi_t4, b_t4 = tend(@phi_s+phi_t3*dt, @b_s+b_t3*dt)
      end
      @phi_s += (phi_t1 + phi_t2*2 + phi_t3*2 + phi_t4 ) * (dt/6.0)
      @b_s += (b_t1 + b_t2*2 + b_t3*2 + b_t4 ) * (dt/6.0)
      if @w_hips
        @hips_s += (hips_t1 + hips_t2*2 + hips_t3*2 + hips_t4 ) * (dt/6.0)
      end
      @step += 1
      @t = dt * @step   # @t + dt にすると丸めが蓄積するので
    end

    def self.n_prof_0(nz)
      npr = NArray.float(nz+1).fill(1.0)
      npr
    end

    # 下端で na, 圏界面以上で nc, そのすぐ下の「TTL 下端」の安定度極小値が nb
    # である N のプロファイルを作る。pow が大きいほど nb になる高度が上がる。
    #
    def self.n_prof_1(nz, nzt, na: 1.2, nb: 0.5, nc: 2.5, pow: 4)
      npr = NArray.float(nz+1)
      x = NArray.float(nzt+1).indgen! * (1.0/nzt)
      d = na-nb
      e = nc-nb
      npr[0..nzt] = d * ( (NMath.sqrt(e/d)+1) * x**pow - 1.0 )**2 + nb
      npr[nzt..-1] = nc if nzt<nz+1   # nztのところでは丸めを除き上式と一致
      npr
    end

    def self.n_prof_2(nz, nzt)
      nltr = 1.0  # 下部対流圏安定度（無次元）：１ぐらいにする
      nst = 2.5   # 成層圏安定度（無次元）
      npr = NArray.float(nz+1)
      npr[0..nzt-1] = nltr
      npr[nzt..-1] = nst
      npr
    end

  end
end

#########################################################
if $0 == __FILE__
  require "numru/ggraph"
  include NumRu
  #nx = 128
  #nx = 256
  #nz = 16
  nx = 1024
  nz = 128

  #nu = kp = 0.0

  nu = kp = 10.0
  nl = 2

  #nu = kp = 10.0
  #nu = kp = 1.0
  #nl = 1

  #linear = true
  linear = false
  #hydro = true
  hydro = false

  #a = 0.1
  a = 0.05

  bm = Boussinesq2D.new(nx,nz, a, n2=1.0, nu, kp, nl,
                        linear: linear, hydro: hydro)
  amp = 5.0
  bm.set_Q(nil, amp)

  #bm.set_Q_4_p(nz-1, 4, fact: 1.05, boostQamp: 1.0, boost_until: 4.0)
  
  #p bm
  phi = bm.phi
  phi[true,1..4] = 2.0
  phi[true,5..-2] = 1.0
  phi[0,1..-2] = 2.0
  x = NArray.float(nx).indgen! * (2*Math::PI/nx)
  phi[true,1] = NMath.sin(x)
  phi[true,2] = NMath.sin(2*x)
  phi[true,3] = NMath.sin(4*x)
  z = (NArray.float(nz-1).indgen!+1) * (Math::PI/(nz-1))
  phi[0,true] = NMath.sin(2*z)
  s = bm.g2s(phi)
  #p s
  phic = bm.s2g(s)
  p (phi-phic).abs.max
  p phi
  phix = bm.s2g_del_x(s)
  p "phi_x", phix
  phiz = bm.s2g_del_z(s)
  p "phi_z", phiz
  p "Q", bm.Q[10..-1,true]
  p "tend", bm.tend
  #dt = 0.05
  dt = 0.1
  DCL.swpset("iwidth",1400)
  DCL.swpset("iheight",600)
  DCL.swpset("lwait",false)
  DCL.swpset("lalt",true)
  DCL.gropn(1)
  DCL.sldiv("y",2,2)
  DCL.sgpset("lfull",true)
  DCL.uzfact(0.7)
  GGraph.set_fig "viewport"=>[0.15,0.85,0.12,0.36]
  GGraph.set_tone "color_bar"=>true
  t = 0
  401.times do |it|
    bm.integ_RK4(dt)
    t = t+dt
    #p "integ", bm.gphi,bm.gb
    if it % 10 == 0
      gphi=bm.gphi
      gb = bm.gb
      GGraph.tone gphi,true,"title"=>gphi.name + "  t=#{sprintf('%g',t)}"
      GGraph.contour gphi,false
      GGraph.tone gb
      GGraph.contour gb,false
      zeta = bm.gzeta
      GGraph.tone zeta
      GGraph.contour zeta,false
      bb = bm.gbb
      GGraph.tone bb, true, "nlev"=>20
      GGraph.contour bb,false, "nlev"=>20
    end
    break if bm.phi.abs.max.to_f > 1e4
    #sleep 0.03
  end
  DCL.grcls
end
