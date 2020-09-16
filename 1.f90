module vars                            !应该除了沙漏部分需要看看其他都定义完成
implicit none
character(len=50),save:: problem
! problem mode
! 0 -> impact
! 1 -> x-ray shock wave

integer(4), save :: problem_mode
integer(4),save:: shape    !0代表平板形状，1代表圆柱壳体形状
integer(4),save:: eos_mode  
real(8),save:: flyer_length
integer(4),save::mat_model(30)   !各层材料类型指示，0->各向同性，1->各向异性

integer(4),save :: yield_mode(30)
! 0 -> rate-unrelated and no strain hardening(Hill_yield)
! 1 -> rate-related and strain hardening 

logical,save::mat0,mat1   !mat0为真时表示各向同性材料，mat1为真时表示各向异性材料

integer(4),save::fail(800000)  !材料单元当前失效判定，0代表材料未失效，1代表材料失效
integer(4),save::fail_his(800000)  !材料单元失效历史，0代表材料未曾失效，1代表材料曾经失效


integer(4),save:: ne,nj
integer(4),save:: nj_xy,nj_xz
integer(4),save::symm_xy(10000),symm_xz(10000)


integer(4),save:: n_cm1,n_cm2,n_cm3
!分别为单元总数，节点总数，碰撞体的节点总数，靶体的节点总数,中线上的节点数
integer(4),save::st,sn  !关于计算步数的控制变量

integer(4),save::num_time_record  !指定输出结果的时刻总数,最多=20
real(8),save::time_record(20)     !指定需要输出结果的具体时刻值

integer(4),save::cm1(8,100000),cm2(8,100000),cm3(1000)  !碰撞问题中，cm1表示碰撞体包含的单元，cm2表示碰撞体包含的单元，cm3表示碰撞接触线所在的(左)单元

integer(4),save::n_cmup  !中轴线上方的单元数量
integer(4),save::cmup(10000)    !沿x方向的中轴上方所包含的单元(平板模型)

integer(4),save::n_cm_r  !任意半径的四分之一圆周上的单元数目
integer(4),save::n_cm_theta !任意角度半径方向上的单元数目
integer(4),save::cm_r1(5000),cm_r2(5000)  !半径为r1和r2的四分之一圆弧上的单元数组
integer(4),save::cm_theta_0(2000),cm_theta_10(2000),cm_theta_20(2000),cm_theta_30(2000),cm_theta_40(2000),cm_theta_45(2000),   &
                 cm_theta_50(2000),cm_theta_60(2000),cm_theta_70(2000),cm_theta_80(2000),cm_theta_90(2000)
!角度分别为0，10，20，30，40，45，50，60，70，80，90的半径上的单元数组

integer(4),parameter::max_layer=30
integer(4),save::lmax   !lmax<max_layer 
integer(4),save::lmax_iso,lmax_aniso     !各向同性材料与各向异性材料的层数  
integer(4),save::ne_layer(800000)  !各单元所属层数

integer(4),save::layer_no_iso(30)   !各向同性材料中第i层材料在整体层数中的位置，主要用于数据输入
integer(4),save::layer_no_aniso(30) !各向异性材料中第i层材料在整体层数中的位置，主要用于数据输入
		
integer(4),save:: ia(9,800000)   !j=1-8表示组成某单元的8个节点编号,j=9代表该单元所在层数的属性
real(8),save:: xyz0(3,800000),xyz(3,800000)  !j=1表示x坐标,j=2表示y坐标,j=3表示z坐标

integer(4),save:: num_elefix_x,num_elefix_y,num_elefix_z  !选取x坐标相同的一组点，选取y坐标相同的一组点，，选取z坐标相同的一组点最多各取20个
integer(4),save:: ne_his_x(20),ne_his_y(20),ne_his_z(20)
integer(4),save:: nscr
!关于输出文件命名的相关变量定义
character(len=18) time_output(20)     !基本形式为time_01_output.txt
character(len=13) fne_his_x(20),fne_his_y(20),fne_his_z(20)         !基本形式为fix_ex_01.dat,fix_ey_01.dat

integer(4) index_ele_write(800000)         !对于圆柱壳体模型，index_ele_write(i)==1代表该单元处于上半平面，需要输出结果

real(8),save:: t,dt,h_dt,t_total   !时间参数
real(8),save:: length_1,length_2,length_3,length_layer(30),dl_1,dl_2,dl_3  
!对于平板模型，分别对应平面x方向总长度，y方向总长度，z方向总长度，各层的厚度，x方向y方向z方向的网格尺寸；
!对于圆柱壳体模型，分别对应圆环的内半径、外半径，各层沿半径方向的厚度，半径方向的网格尺寸，角度划分尺寸
real(8),save::length_x,length_y,length_z,lx(30),dx,dy,dz
!对应平面x方向总长度，y方向总长度，z方向总长度各层的厚度，x方向y方向z方向的网格尺寸
real(8),save:: r_in,r_out,r_height,lr(30),dr,d_theta,d_height     
!分别为圆环的内半径、外半径，各层沿半径方向的厚度，半径方向的网格尺寸,角度划分尺寸



real(8),save::rho0(30)  !各层的初始密度
real(8),save::rho0_iso(30),rho0_aniso(30)

real(8),save::c11_0(30),c22_0(30),c33_0(30),c12_0(30),c13_0(30),c23_0(30),c44_0(30),c55_0(30),c66_0(30)         !初始各层本构系数，由刚度矩阵转换得到,包括3*3的三个正应力关系和3个剪切关系
real(8),save::c123_0(6,6,20)
real(8),save::c11(800000),c22(800000),c33(800000),c12(800000),c21(800000),c13(800000),c31(800000),c23(800000),c32(800000),c44(800000),c55(800000),c66(800000)         !各单元本构系数
real(8),save::c123(6,6,800000)

real(8),save::ex_0(30),ey_0(30),ez_0(30),gxy_0(30),gxz_0(30),gyz_0(30)!,pxy(30),pxz(30),pyz(30)!三个弹性模量，剪切模量，泊松比
real(8),save::ex_aniso(30),ey_aniso(30),ez_aniso(30),pxy_aniso(30),pxz_aniso(30),pyz_aniso(30),gxy_aniso(30),gxz_aniso(30),gyz_aniso(30)
real(8),save::ex_iso(30),ey_iso(30),ez_iso(30),pxy_iso(30),pxz_iso(30),pyz_iso(30),gxy_iso(30),gxz_iso(30),gyz_iso(30)
real(8),save::ex(800000),ey(800000),ez(800000),gxy(800000),gxz(800000),gyz(800000),pxy(800000),pxz(800000),pyz(800000)

real(8),save::cos_theta(800000),sin_theta(800000)  !圆柱壳体模型中，初始各单元材料主轴(各向异性材料)与坐标轴所成夹角的余弦值和正弦值

real(8),save::strs_yld_0(30),strs_yld(800000)!各向同性材料的屈服强度,前者为以层表示，后者以单元表示

real(8),save::strs_yld_110(30),strs_yld_220(30),strs_yld_330(30),strs_yld_120(30),strs_yld_130(30),strs_yld_230(30)  
real(8),save::strs_yld_11(800000),strs_yld_22(800000),strs_yld_33(800000),strs_yld_12(800000),strs_yld_13(800000),strs_yld_23(800000)
real(8),save::strs_yld_123(3,3,800000)
!各向异性材料的屈服强度,前者为以层表示，后者以单元表示
real(8),save::f1_yield(800000)


real(8),save::strs_frct_0(30),strs_frct(800000)
real(8),save::strs_frct_123(3,3,800000)

real(8),save::strs_frct_110(30),strs_frct_220(30),strs_frct_330(30),strs_frct_120(30),strs_frct_130(30),strs_frct_230(30)
real(8),save::strs_frct_11(800000),strs_frct_22(800000),strs_frct_33(800000),strs_frct_12(800000),strs_frct_13(800000),strs_frct_23(800000)

integer(4),save::gas(800000)  !0代表物质不是气体，1代表物质为气体

real(8),save::c0(30),s(30),gama(30)  !物态方程参数
real(8),save::c0_iso(30),s_iso(30),gama_iso(30)
real(8),save::c0_aniso(30),s_aniso(30),gama_aniso(30)

real(8),save::h0(30)      !PUFF物态方程参数
real(8),save::h0_iso(30) 
real(8),save::h0_aniso(30) 

real(8),save::cn(30)

real(8),save::es(30)         !升华能
real(8),save::es_iso(30)
real(8),save::es_aniso(30)

real(8),save::q1,q2    !粘性系数
real(8),save::qhg   !沙漏粘性系数 取为0.1

real(8),save::vf

real(8),save::ax(30),ay(30),az(30),nx(30),ny(30),nz(30),axy(30),nxy(30),axz(30),nxz(30),ayz(30),nyz(30)
real(8),save::ax_aniso(30),ay_aniso(30),az_aniso(30),nx_aniso(30),ny_aniso(30),nz_aniso(30),axy_aniso(30),nxy_aniso(30),axz_aniso(30),nxz_aniso(30),ayz_aniso(30),nyz_aniso(30)
real(8),save::ax_iso(30),nx_iso(30)

real(8),save::rate(800000),strain_rate_eff(800000)
real(8),save::beta0(30),strain_rate0(30)        !应变率有效因子的参数0.02183,初始参考应变率 取值为0.001/s,应变率因子
real(8),save::beta0_iso(30),strain_rate0_iso(30)
real(8),save::beta0_aniso(30),strain_rate0_aniso(30)

integer(4),save::plastic_index(800000)   !塑性变形指示因子，为0时未发生塑性变形，为1时发生塑性变形


real(8),save::u(3,800000),du(3,800000),v(3,800000),a(3,800000),v_mid(3,800000)  !在3nj个自由度上的节点位移,节点速度,节点加速度,中间时刻节点速度
real(8),save::mass(800000)     !在nj个自由度上的质量矩阵（对角化）
real(8),save::x0(800000),y0(800000),z0(800000)   !初始时每单元中心的整体坐标
real(8),save::x(800000),y(800000),z(800000),v_cen_x(800000),v_cen_y(800000),v_cen_z(800000),l,h,m      !每单元中心的整体坐标，速度,自然坐标
!real(8),save::v_cen_mid_x(800000),v_cen_mid_y(800000)    !每单元中心处在整体坐标中中间时刻的速度
!real(8),save::v_cen_mid_1(800000),v_cen_mid_2(800000)    !每单元中心处在主轴坐标中中间时刻的速度
! real(8),save::velocity_xyz(3,8,800000)
real(8),save::shape_n(8),shape_n0(8) !形状函数,中心点的形状函数
real(8),save::dn_l(8),dn_h(8),dn_m(8),dn0_l(8),dn0_h(8),dn0_m(8)  !形状函数、中心点的形状函数对自然坐标l,h,m的导数
real(8),save::dn_lhm(8,3)
real(8),save::dn_xyz(8,3)
real(8),save::coord_xyz(3,8)
real(8),save::coord_xyz0(3,8)
real(8),save::dn_x(8),dn_y(8),dn_z(8)           !形状函数对整体坐标x,y,z的导数
real(8),save::jj(3,3),jj_inverse(3,3)       !(x,y,z)对(l,h,m)偏导数的雅各比矩阵,矩阵行列式，逆矩阵
real(8),save::jj0(3,3),jj0_inverse(3,3)    !(x0,y0,z0)对(l,h,m)偏导数的雅各比矩阵,矩阵行列式，逆矩阵
real(8),save::f11,f12,f13,f21,f22,f23,f31,f32,f33                  !(x,y,z)对(x0,y0,z0)的偏导数的雅各比矩阵，变形梯度矩阵
real(8),save::deform(3,3)
real(8),save::rotate_R(3,3)
real(8),save::rotate_original(3,3)


real(8),save::rho(800000),rho_old(800000),volume0(800000),volume(800000),volume_old(800000),dvolume(800000),volume_average(800000)
real(8),save::vol0(800000),vol(800000),vol_old(800000),dvol(800000),vol_average(800000)  !单元密度，现时面积与上一时刻的面积
real(8),save::e(800000),de(800000),p(800000),dp(800000),ee         !各单元的内能,平均应力,总内能

real(8),save::np(800000)

real(8),save::p_dst,p_dsp         !偏应变贡献压力，塑性应变贡献压力

real(8),save::f_ext(3,800000),f_int(3,800000),f_ext_old(3,800000),f_int_old(3,800000) 


real(8),save::stress_xx(800000),stress_yy(800000),stress_zz(800000),stress_xy(800000),stress_xz(800000),stress_yz(800000) !每个单元（l=0,s=0,m=0）处的应力
real(8),save::stress_xyz(3,3,80000)
real(8),save::stress_11(800000),stress_22(800000),stress_33(800000),stress_12(800000),stress_13(800000),stress_23(800000) !在主轴坐标中，每个单元（l=0,s=0,m=0）处的应力
real(8),save::stress_123(3,3,80000)
real(8),save::stress_d_11(800000),stress_d_22(800000),stress_d_33(800000),stress_d_12(800000),stress_d_13(800000),stress_d_23(800000)   !每个单元（l=0,s=0,m=0）处当前时刻的偏应力
real(8),save::stress_d_123(3,3,800000)

real(8),save::stress_node_xx(800000),stress_node_yy(800000),stress_node_zz(800000),stress_node_xy(800000),stress_node_xz(800000),stress_node_yz(800000),p_node(800000),e_node(800000) 
real(8),save::stress_node_xyz(3,3,80000)

real(8),save::diag_vector(3)
real(8),save::diag_matrix(3,3)
real(8),save::diag_sum
integer,save::number1,number2,number3
!球—偏分解专用变量

real(8),save::q(800000),q_old(800000),q_av(800000)   !人工粘性项（当前时刻和上一时刻的值及平均值）
!==========================================================================================三维的沙漏修正
real(8),save::hggama1(800000),hggama2(800000),hggama3(800000),hggama4(800000)  !每个单元沙漏基矢量γ的四个分量（三维情况下四个是否够？？）
real(8),save::f_hg(3,800000)     !每个节点处的x,y方向的沙漏力

real(8),save::strain_xx(800000),strain_yy(800000),strain_zz(800000),strain_xy(800000),strain_xz(800000),strain_yz(800000)      !每个单元（l=0,s=0,m=0）处的应变,等效塑性应变率
real(8),save::strain_xyz(3,3,80000)

real(8),save::strain_xx_rate(800000),strain_yy_rate(800000),strain_zz_rate(800000),strain_xy_rate(800000),strain_xz_rate(800000),strain_yz_rate(800000)   !应变率分量
real(8),save::strain_xyz_rate(3,3,800000)
real(8),save::strain_11(800000),strain_22(800000),strain_33(800000),strain_12(800000),strain_13(800000),strain_23(800000)      !在主轴坐标中，每个单元（l=0,s=0,m=0）处的应变
real(8),save::strain_123(3,3,800000)

real(8),save::strain_11_rate(800000),strain_22_rate(800000),strain_33_rate(800000),strain_12_rate(800000),strain_13_rate(800000),strain_23_rate(800000)   !主轴坐标系中应变率分量
!real(8),save::strain_123_rate(3,3,800000)
real(8),save::strain_p_11(800000),strain_p_22(800000),strain_p_33(800000),strain_p_12(800000),strain_p_13(800000),strain_p_23(800000)      !在主轴坐标中，每个单元（l=0,s=0,m=0）处的塑性应变
real(8),save::strain_p_123(3,3,800000)

real(8),save::strain_p_eff(800000),strain_p_rate_eff(800000)  !等效塑性应变,等效塑性应变增量，等效塑性应变率

!real(8),save::peak_axis(800000)  !中线上各单元的应力峰值

real(8),save::rad_energy_total  !入射通量密度
integer(4),save::rad_time_mode  !时间谱模式，0为矩形谱，1为三角形谱
real(8),save::rad_time_width    !时间谱的宽度，三角形谱为半高宽
integer(4),save::rad_energy_num !时间份数=时间谱宽度/dt
real(8),save::rad_energy(8000000)   !功率密度(energy_num维)

integer(4),parameter::ken=5000
real(8),save::wl_short,wl_long
real(8),save::wl(5000)
real(8),save::wer(5000)

integer(4),save::n_blackbody  !黑体谱的个数，最多取20组
real(8),save::t_kev(20)       !各黑体组分的黑体温度，单位为keV
real(8),save::t_k(20)       !各黑体组分的黑体温度，单位为k(开)
real(8),save::b_energy(20)    !各黑体组分所占总能量的能量份额
real(8),save::f_blackbody_each(20)  !各黑体组分的普朗克函数
real(8),save::a_planck(20)    !各黑体谱辐射强度的比例系数

real(8),save::energy_der(800000)
real(8),save::energy_xray_total(800000)
real(8),save::energy_xray(800000)

integer(4),parameter:: max_element_num=10

real(8),save::parameter_A1,parameter_A2,parameter_A3,parameter_B0,parameter_T1                  
real(8),save::a_11,a_22,a_33,a_12,a_13,a_23,a_44,a_55,a_66,sigma_well_1,sigma_well_10,mean_plastic    
real(8),save::tensil_11,tensil_22,tensil_33,tensil_12,tensil_23,tensil_31                        
real(8),save::fracture_e11,fracture_e22,fracture_e33,fracture_e12,fracture_e23,fracture_e31     
real(8),save::damage_coupling


type::layer_xray
    integer(4) element_num    !各层的元素总数，<=10

	integer(4) element_z(max_element_num)    !元素的原子序数，元素种类不超过10个
	real(8) element_wp(max_element_num)     !各元素的质量百分比

	real(8) mu(ken)
end type layer_xray

type(layer_xray), save :: xlayer(max_layer)

real(8),save::w_int,w_ext,w_kin,w_orgn,w_err     !内能，外能，动能，初始能量，能量误差百分比
!体积计算相关变量
real,save::core(3,1)
real,save::volume_s(6)






end module vars   