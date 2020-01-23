# -*- coding: utf-8 -*-
=begin
= 無次元化した Anelastic 鉛直ー水平２次元モデル
x: cyclic, 範囲 0..2pi
z: bounded, 範囲 0..pi

Boussinesq2D を継承して定義

== 支配方程式
  η = exp(c*z) (c=Z/(2H)
  ζ_t + αη(u ζ_x+w (ζ_z+2cζ) = -αb_x - ν((-Δ)^nl)ζ
  b_t + αη(u b_x+w b_z) + αN2 w = αQ - κ((-Δ)^nl)ζ
  ζ = α^2 φ_xx + φ_zz + c**2 φ,  u=φ_z-cφ, w=-φ_x


=end

require File.dirname(File.expand_path(__FILE__)) + "/boussinesq2D.rb"

module NumRu
  class Anelastic2D < Boussinesq2D
    def initialize(*args)
      # 追加キーワード引数:
      @nondim_H = args[-1].delete(:nondim_H)   # H/Z

      super(*args)
      raise("need to specify nondim_H (non dim scale hight)") unless @nondim_H
      @cs = 0.5/@nondim_H  # Z/(2*H)
      
      @lapla = @lapla - @cs**2
      @ilapla = 1.0/@lapla

      @eta = NMath.exp(@cs*@ze).newdim(0)
    end

    def gphi
      @gphi = super
      @gphi *= @eta if @dimensional_out
      @gphi
    end
    def gb
      @gb = super
      @gb *= @eta if @dimensional_out
      @gb
    end
    %w!bb hips zeta u w q!.each do |vn|
      eval <<-EOS
        def g#{vn}
          g = super
          g *= @eta if @dimensional_out
          g
        end
      EOS
    end

    def phi_s2u_e(phi_s)
      u = s2g_del_z(phi_s) - g_ze(@cs * phi_s)
    end
    def phi_s2u(phi_s)
      u = g_zc( s2g_del_z(phi_s) ) - @cs * phi_s
    end
    def zeta_z_s_4adv(zeta_s)
      g_zc(s2g_del_z(zeta_s)) + (2*@cs) * s2g(zeta_s)
    end
    def adv_factor
      @a*@eta[true,1..-2]
    end

  end
end
