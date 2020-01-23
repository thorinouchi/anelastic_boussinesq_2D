# -*- coding: utf-8 -*-
=begin
= 対流バーストに伴う重力波の簡単モデル計算
(2D, x-z座標)


== 実行例

nohup ruby boussinesq2D_CBmodel.rb --dt 0.05 --skip 6 --zrange 32 --n0 0.01 --alpha 0.08 --amp 5 --nl 2 --nu 10 --kp 10 --nx 512 --nprof 2 --nQ 2 --sponge_nl 6 --sponge_rate 5 --tend 72.1


=end


require File.dirname(File.expand_path(__FILE__)) + "/boussinesq2D.rb"
require File.dirname(File.expand_path(__FILE__)) + "/anelastic2D.rb"
require "numru/ggraph"
require "optparse"
include NumRu

command = "ruby #{$0} "+ ARGV.join(" ")

opt = OptionParser.new

modeltype = Boussinesq2D
opt.on("--anelastic"){ modeltype = Anelastic2D}

dt = 0.1  # dt for time integration
opt.on("--dt VAL",Float){|x| dt=x}

tend = 30
opt.on("--tend VAL",Float){|x| tend=x}

opath = nil  # ファイル名を指定したら出力
opt.on("-o PATH"){|x| opath=x}

fname = nil  # 画像ファイル名prefix
opt.on("--fname PREFIX"){|x| DCL.swpset("fname",x)}

skip = 5  # -> 出力時間間隔 ＝ skip * dt
opt.on("--skip VAL",Integer){|x| skip=x}

nx = 256   # 水平格子点数
opt.on("--nx VAL",Integer){|x| nx=x}

nz = 64    # 鉛直層数
opt.on("--nz VAL",Integer){|x| nz=x}

nzt = nil   # 対流圏の層数 (default: nz/2)
opt.on("--nzt VAL",Integer){|x| nzt=x}

n_ul = nil   #  加熱の上限レベル (nz-1以下)
opt.on("--n_ul VAL",Integer){|x| n_ul=x}

nprof = 1   # N profile の番号 (0, 1, 2)
opt.on("--nprof VAL",Integer){|x| nprof=x}

sponge_nl = nil  # 整数ならスポンジあり
opt.on("--sponge_nl VAL",Integer){|x| sponge_nl=x}

sponge_rate = 1.0  # 最大(最上層で)の時定数
opt.on("--sponge_rate VAL", Float){|x| sponge_rate=x}

nQ = 3   # heating番号 (2, 3, 4)
opt.on("--nQ VAL",Integer){|x| nQ=x}

amp = 1.0  # heating amplitude (when nQ = 2,3)
opt.on("--amp VAL",Float){|x| amp=x}

amp2 = nil  # wave-2 amplitude (when nQ = 2)
opt.on("--amp2 VAL",Float){|x| amp2=x}

pq_fact = 1.05  # for nQ=4 (簡易積雲P): 条件付き不安定度ファクター
opt.on("--pq_fact VAL",Float){|x| pq_fact = x}

pq_b_amp = 0.5  # for nQ=4 (簡易積雲P): 初期boost用Qの最大値
opt.on("--pq_b_amp VAL",Float){|x| pq_b_amp = x}

pq_b_until = 1  # for nQ=4 (簡易積雲P): 初期boost適用いつまで
opt.on("--pq_b_until VAL",Float){|x| pq_b_until = x}

nw = 6   # heating width
opt.on("--nw VAL",Integer){|x| nw=x}

nl = 2  # 数値粘性/拡散の次数
opt.on("--nl VAL",Integer){|x| nl=x}

nu = 10.0  # 数値粘性の係数 (最大波数の減衰率)
opt.on("--nu VAL",Float){|x| nu=x}

kp = 10.0  # 数値拡散の係数 (最大波数の減衰率)
opt.on("--kp VAL",Float){|x| kp=x}

k_v_h = 1.0  # 鉛直拡散粘性を k_v_h に調整
opt.on("--k_v_h VAL",Float){|x| k_v_h=x}

alpha = 0.1  # 縦/横
opt.on("--alpha VAL",Float){|x| alpha=x}

zrange = nil  # [km] 例 10. 高度範囲（次元付き読み替え用：計算には影響なし）
opt.on("--zrange VAL",Float){|x| zrange=x}

scale_height = 8.0  # [km] Anelastic な場合のスケールハイト 
opt.on("--scale_height VAL",Float){|x| scale_height=x}

n0 = nil  # 例 0.01. 基準浮力振動数（次元付き読み替え用：計算には影響なし）
opt.on("--n0 VAL",Float){|x| n0=x}

linear = false  # 線形（移流項無視）
opt.on("--linear"){linear = true}

hydro = false  # 静水圧強制（linearのときのみ意味がある）
opt.on("--hydro"){hydro = true}

w_hips = true  # heat-induced passive scalar あり
opt.on("--wo_hips"){w_hips = false}

iws = 1
opt.on("--png"){
  DCL.swpset('ifl',1)
  iws=2
}

opt.parse!(ARGV)

nzt = nz/2 unless nzt

case nprof
when 0
  n = modeltype.n_prof_0(nz)
  n_ul = nz unless n_ul
when 1
  np1_na = 1.2
  np1_nb = 0.5
  np1_nc = 2.5
  np1_pow = 4
  n = modeltype.n_prof_1(nz, nzt, na: np1_na, nb: np1_nb, nc: np1_nc, pow: np1_pow)
  unless n_ul
    # 加熱の上限をnが極小となる層までにする
    k = nzt
    while( n[k]-n[k-1] > 0.0)
      k -= 1
      break if k<0
    end
    n_ul = k
  end
when 2
  n = modeltype.n_prof_2(nz, nzt)
  n_ul = nzt unless n_ul    # 加熱の上限 (nzt以下であるべき)
else
  raise("undfined N-profile number: #{nprof}")
end

#puts "** N"
#(0...n.length).each{|i| printf "%d  %f  %f\n", i, 2.0*(i+1)/nz,n[i]}

n2 = n.newdim(0)**2

hs = zrange ? zrange*1000/Math::PI : nil
ts = n0 ? 1.0/n0 : nil

nondim_H = (scale_height*1000) / hs

bm = modeltype.new(nx,nz, alpha, n2, nu, kp, nl,
                      linear: linear, hydro: hydro, w_hips: w_hips,
                      hs: hs, ts: ts, k_v_h: k_v_h, nondim_H: nondim_H)
case nQ
when 2
  bm.set_Q_2(n_ul, nw, amp, amp2)
when 3
  bm.set_Q_3(n_ul, nw, amp)
when 4
  bm.set_Q_4_p(n_ul, nw, fact: pq_fact, boostQamp: pq_b_amp,
               boost_until: pq_b_until)
else
  raise("unsupported heating number: #{nQ}")
end

if sponge_nl
  bm.set_sponge(sponge_nl, sponge_rate)
end

DCL.swpset("iwidth",1600)
DCL.swpset("iheight",800)
DCL.swpset("lwait",false)
DCL.swpset("lalt",true)
DCL.gropn(iws)
DCL.sgscmn(4)
DCL.sldiv("y",2,3)
DCL.sgpset("lfull",true)
DCL.uzfact(0.55)
GGraph.set_fig "viewport"=>[0.1,0.9,0.12,0.27]
nitr = (tend/dt).round
#gint = [(0.2/dt).round,1].max

if opath
  ts = Array.new
  phis = Array.new
  bs = Array.new
  bbs = Array.new
  zetas = Array.new
  hipss = Array.new
  us = Array.new
  ws = Array.new
  qs = Array.new
end

cbvlen = 0.14

(0...nitr).each do |istp|
  #puts "#{istp+1} / #{nitr}"
  if (istp % skip == 0)
    gphi=bm.gphi
    GGraph.tone gphi,true,"title"=>gphi.name + "  t=#{sprintf('%g',bm.t_d)}"
    GGraph.contour gphi,false
    GGraph.color_bar "vlen"=>cbvlen

    gb = bm.gb
    GGraph.tone gb
    GGraph.contour gb,false
    GGraph.color_bar "vlen"=>cbvlen

    zeta = bm.gzeta
    GGraph.tone zeta
    #GGraph.contour zeta,false
    GGraph.color_bar "vlen"=>cbvlen

    hips = bm.ghips
    if hips.max==0
      GGraph.tone hips, true, "clr_min"=>54, "min"=>0,
                "max"=>amp, "inf_max"=>true
      #GGraph.color_bar "vlen"=>cbvlen
    else
      if amp > 1.0
        min = amp*1e-3
      else
        min = 1e-3
      end
      GGraph.tone hips, true, "clr_min"=>54,
                "log"=>true, "min"=>min
      GGraph.color_bar "vlen"=>cbvlen, "log"=>true
    end

    bb = bm.gbb
    GGraph.tone bb, true, "nlev"=>20
    if modeltype==Boussinesq2D
      max = 3
    else
      max = 9
    end
    GGraph.contour bb,false, "nlev"=>15, "max"=>max
    GGraph.color_bar "vlen"=>cbvlen

    gu = bm.gu
    GGraph.tone gu, true, "title"=>"u & w"
    gw = bm.gw
    GGraph.contour gw, false
    GGraph.color_bar "vlen"=>cbvlen
    
    if opath
      ts.push(bm.t_d)
      phis.push(gphi.copy)
      bs.push(gb.copy)
      bbs.push(bb)
      zetas.push(zeta)
      hipss.push(hips.copy)
      us.push(gu)
      ws.push(gw)
      gq = bm.gq
      qs.push(gq)
    end

  end
  #sleep 0.03

  bm.integ_RK4(dt)
  if bm.phi.abs.max.to_f > 1e4
    puts "** It's likely that instability has occured. Stopping."
    break 
  end
end

#loop do
#  DCL.swflsh
#  sleep 0.5
#end

if opath
  tattr = {"long_name"=>"time","units"=>bm.t_d_un}
  phi = GPhys.concat(phis, ts, "t", tattr)
  b = GPhys.concat(bs, ts, "t", tattr)
  bb = GPhys.concat(bbs, ts, "t", tattr)
  hips = GPhys.concat(hipss, ts, "t", tattr)
  zeta = GPhys.concat(zetas, ts, "t", tattr)
  u = GPhys.concat(us, ts, "t", tattr)
  w = GPhys.concat(ws, ts, "t", tattr)
  q = GPhys.concat(qs, ts, "t", tattr)

  ofl = NetCDF.create(opath)
  ofl.put_att("command", command)
  ofl.put_att("history", "#{Time.now}; #{ENV['USER']}; #{__FILE__}")
  v_alpha = ofl.def_var("alpha", "float", [])
  v_alpha.put_att("long_name","aspect ratio (dz/dx)")
  v_alpha.put_att("units","")
  if hs  # dimensional out put
    v_zrange = ofl.def_var("zrange", "float", [])
    v_zrange.put_att("long_name","height range")
    v_zrange.put_att("units","km")
    v_N0 = ofl.def_var("N0", "float", [])
    v_N0.put_att("long_name","reference N")
    v_N0.put_att("units","s-1")
    ofl.enddef
    v_alpha.put(alpha)
    v_zrange.put(zrange)
    v_N0.put(n0)
    ofl.redef
    GPhys::IO.write_grid(ofl,phi.grid)
    v_n = ofl.def_var("n", "float", ["z"])
    v_n.put_att("long_name","buoyancy frequency")
    v_n.put_att("units","s-1")
    ofl.enddef
    v_n.put(NMath.sqrt(bm.n2e)*n0)
  else
    ofl.enddef
    v_alpha.put(alpha)
  end
  ofl.redef
  GPhys::IO.write(ofl,phi)
  GPhys::IO.write(ofl,b)
  GPhys::IO.write(ofl,hips)
  GPhys::IO.write(ofl,bb)
  GPhys::IO.write(ofl,zeta)
  GPhys::IO.write(ofl,u)
  GPhys::IO.write(ofl,w)
  GPhys::IO.write(ofl,q)
  ofl.close
end

DCL.grcls
