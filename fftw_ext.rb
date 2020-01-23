# -*- coding: utf-8 -*-
require "narray"
require "numru/fftw3"

module NumRu
  module FFTW3
    @@fft_forward = -1
    @@fft_backward = 1

    module_function

    # sine 変換 (ただし出力は実部がゼロの複素 NArray)．特定の次元についてのみ．
    # 
    # 引数
    # * g [NArray] 実NArray (第dim次元では両端のゼロを除く)
    # * dim [Integer] 変換する次元
    # 
    # 出力
    # * f [NArray]
    def sin_fw(g, dim=0)
      raise(ArgumentError, "Invalid dim #{dim} >= #{g.rank}") if dim >= g.rank
      n = 2*(g.shape[dim]+1)  # 正規化係数は複素の場合の最大波数
      f = FFTW3.fft_r2r(g, FFTW3::RODFT00, dim) * (1.0/n)
      f
    end

    # sin_fw の逆変換 (入力は実部がゼロの複素 NArray, 出力は実NArray)
    def sin_bk(f, dim=0)
      raise(ArgumentError, "Invalid dim #{dim} >= #{f.rank}") if dim >= f.rank
      FFTW3.fft_r2r(f, FFTW3::RODFT00, dim)
    end

    # sin_fw の出力(両端がない) をcosの係数とみなして逆変換 （微分用）
    # (定数項と最大波数項はゼロとする)
    # dim次元は2個長くなる。
    def cos_bk_sinout(f, dim=0)
      raise(ArgumentError, "Invalid dim #{dim} >= #{f.rank}") if dim >= f.rank
      sh = f.shape
      n = sh[dim]
      sh[dim] = n+2
      fe = NArray.float(*sh)
      fe[ *( [true]*dim + [1..n,false] ) ] = f
      FFTW3.fft_r2r(fe, FFTW3::REDFT00, dim)
    end

    # 周期境界関数の複素 FFT forward
    # （入力は実/複素どちらでもいいが, 主に実数を想定）
    def fft_fw(g, dim=0)
      n = g.shape[dim]
      FFTW3.fft_r2r(g, FFTW3::R2HC, dim) * (1.0/n)
    end

    # 周期境界関数の複素 FFT backward
    # （もとの入力が実数であることを想定し，確認せず実部のみとって返す）
    def fft_bk(f, dim=0)
      FFTW3.fft_r2r(f, FFTW3::HC2R, dim)
    end

  end
end

########################################
if $0 == __FILE__
  include NumRu

  eps = 1e-13

  zg = NArray.float(10,6,2)
  zg[true,0,true] = 1
  zg[true,-1,true] = -1

  zg.indgen!
  zg[true,0,true] = 0
  zg[true,-1,true] = 0
  puts "\n** test sine transform"
  p zg
  zf = FFTW3.sin_fw(zg,1)
  p zf
  zgb = FFTW3.sin_bk(zf,1)
  difmax = (zgb - zg).abs.max
  if difmax < eps
    print "Test passed:  max residual= #{difmax}\n"
  else
    raise "Test failed:  max residual= #{difmax}"
  end

  zgc = FFTW3.cos_bk_sinout(zf,1)
  p zgc
  
  puts "\n** test preiodic FFT"
  zf = FFTW3.fft_fw(zg,1)
  zgb = FFTW3.fft_bk(zf,1)
  difmax = (zgb - zg).abs.max
  if difmax < eps
    print "Test passed:  max residual= #{difmax}\n"
  else
    raise "Test failed:  max residual= #{difmax}"
  end

end

