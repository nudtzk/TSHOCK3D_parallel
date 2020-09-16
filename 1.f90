module vars                            !Ӧ�ó���ɳ©������Ҫ�����������������
implicit none
character(len=50),save:: problem
! problem mode
! 0 -> impact
! 1 -> x-ray shock wave

integer(4), save :: problem_mode
integer(4),save:: shape    !0����ƽ����״��1����Բ��������״
integer(4),save:: eos_mode  
real(8),save:: flyer_length
integer(4),save::mat_model(30)   !�����������ָʾ��0->����ͬ�ԣ�1->��������

integer(4),save :: yield_mode(30)
! 0 -> rate-unrelated and no strain hardening(Hill_yield)
! 1 -> rate-related and strain hardening 

logical,save::mat0,mat1   !mat0Ϊ��ʱ��ʾ����ͬ�Բ��ϣ�mat1Ϊ��ʱ��ʾ�������Բ���

integer(4),save::fail(800000)  !���ϵ�Ԫ��ǰʧЧ�ж���0�������δʧЧ��1�������ʧЧ
integer(4),save::fail_his(800000)  !���ϵ�ԪʧЧ��ʷ��0�������δ��ʧЧ��1�����������ʧЧ


integer(4),save:: ne,nj
integer(4),save:: nj_xy,nj_xz
integer(4),save::symm_xy(10000),symm_xz(10000)


integer(4),save:: n_cm1,n_cm2,n_cm3
!�ֱ�Ϊ��Ԫ�������ڵ���������ײ��Ľڵ�����������Ľڵ�����,�����ϵĽڵ���
integer(4),save::st,sn  !���ڼ��㲽���Ŀ��Ʊ���

integer(4),save::num_time_record  !ָ����������ʱ������,���=20
real(8),save::time_record(20)     !ָ����Ҫ�������ľ���ʱ��ֵ

integer(4),save::cm1(8,100000),cm2(8,100000),cm3(1000)  !��ײ�����У�cm1��ʾ��ײ������ĵ�Ԫ��cm2��ʾ��ײ������ĵ�Ԫ��cm3��ʾ��ײ�Ӵ������ڵ�(��)��Ԫ

integer(4),save::n_cmup  !�������Ϸ��ĵ�Ԫ����
integer(4),save::cmup(10000)    !��x����������Ϸ��������ĵ�Ԫ(ƽ��ģ��)

integer(4),save::n_cm_r  !����뾶���ķ�֮һԲ���ϵĵ�Ԫ��Ŀ
integer(4),save::n_cm_theta !����ǶȰ뾶�����ϵĵ�Ԫ��Ŀ
integer(4),save::cm_r1(5000),cm_r2(5000)  !�뾶Ϊr1��r2���ķ�֮һԲ���ϵĵ�Ԫ����
integer(4),save::cm_theta_0(2000),cm_theta_10(2000),cm_theta_20(2000),cm_theta_30(2000),cm_theta_40(2000),cm_theta_45(2000),   &
                 cm_theta_50(2000),cm_theta_60(2000),cm_theta_70(2000),cm_theta_80(2000),cm_theta_90(2000)
!�Ƕȷֱ�Ϊ0��10��20��30��40��45��50��60��70��80��90�İ뾶�ϵĵ�Ԫ����

integer(4),parameter::max_layer=30
integer(4),save::lmax   !lmax<max_layer 
integer(4),save::lmax_iso,lmax_aniso     !����ͬ�Բ�����������Բ��ϵĲ���  
integer(4),save::ne_layer(800000)  !����Ԫ��������

integer(4),save::layer_no_iso(30)   !����ͬ�Բ����е�i���������������е�λ�ã���Ҫ������������
integer(4),save::layer_no_aniso(30) !�������Բ����е�i���������������е�λ�ã���Ҫ������������
		
integer(4),save:: ia(9,800000)   !j=1-8��ʾ���ĳ��Ԫ��8���ڵ���,j=9����õ�Ԫ���ڲ���������
real(8),save:: xyz0(3,800000),xyz(3,800000)  !j=1��ʾx����,j=2��ʾy����,j=3��ʾz����

integer(4),save:: num_elefix_x,num_elefix_y,num_elefix_z  !ѡȡx������ͬ��һ��㣬ѡȡy������ͬ��һ��㣬��ѡȡz������ͬ��һ�������ȡ20��
integer(4),save:: ne_his_x(20),ne_his_y(20),ne_his_z(20)
integer(4),save:: nscr
!��������ļ���������ر�������
character(len=18) time_output(20)     !������ʽΪtime_01_output.txt
character(len=13) fne_his_x(20),fne_his_y(20),fne_his_z(20)         !������ʽΪfix_ex_01.dat,fix_ey_01.dat

integer(4) index_ele_write(800000)         !����Բ������ģ�ͣ�index_ele_write(i)==1����õ�Ԫ�����ϰ�ƽ�棬��Ҫ������

real(8),save:: t,dt,h_dt,t_total   !ʱ�����
real(8),save:: length_1,length_2,length_3,length_layer(30),dl_1,dl_2,dl_3  
!����ƽ��ģ�ͣ��ֱ��Ӧƽ��x�����ܳ��ȣ�y�����ܳ��ȣ�z�����ܳ��ȣ�����ĺ�ȣ�x����y����z���������ߴ磻
!����Բ������ģ�ͣ��ֱ��ӦԲ�����ڰ뾶����뾶�������ذ뾶����ĺ�ȣ��뾶���������ߴ磬�ǶȻ��ֳߴ�
real(8),save::length_x,length_y,length_z,lx(30),dx,dy,dz
!��Ӧƽ��x�����ܳ��ȣ�y�����ܳ��ȣ�z�����ܳ��ȸ���ĺ�ȣ�x����y����z���������ߴ�
real(8),save:: r_in,r_out,r_height,lr(30),dr,d_theta,d_height     
!�ֱ�ΪԲ�����ڰ뾶����뾶�������ذ뾶����ĺ�ȣ��뾶���������ߴ�,�ǶȻ��ֳߴ�



real(8),save::rho0(30)  !����ĳ�ʼ�ܶ�
real(8),save::rho0_iso(30),rho0_aniso(30)

real(8),save::c11_0(30),c22_0(30),c33_0(30),c12_0(30),c13_0(30),c23_0(30),c44_0(30),c55_0(30),c66_0(30)         !��ʼ���㱾��ϵ�����ɸնȾ���ת���õ�,����3*3��������Ӧ����ϵ��3�����й�ϵ
real(8),save::c123_0(6,6,20)
real(8),save::c11(800000),c22(800000),c33(800000),c12(800000),c21(800000),c13(800000),c31(800000),c23(800000),c32(800000),c44(800000),c55(800000),c66(800000)         !����Ԫ����ϵ��
real(8),save::c123(6,6,800000)

real(8),save::ex_0(30),ey_0(30),ez_0(30),gxy_0(30),gxz_0(30),gyz_0(30)!,pxy(30),pxz(30),pyz(30)!��������ģ��������ģ�������ɱ�
real(8),save::ex_aniso(30),ey_aniso(30),ez_aniso(30),pxy_aniso(30),pxz_aniso(30),pyz_aniso(30),gxy_aniso(30),gxz_aniso(30),gyz_aniso(30)
real(8),save::ex_iso(30),ey_iso(30),ez_iso(30),pxy_iso(30),pxz_iso(30),pyz_iso(30),gxy_iso(30),gxz_iso(30),gyz_iso(30)
real(8),save::ex(800000),ey(800000),ez(800000),gxy(800000),gxz(800000),gyz(800000),pxy(800000),pxz(800000),pyz(800000)

real(8),save::cos_theta(800000),sin_theta(800000)  !Բ������ģ���У���ʼ����Ԫ��������(�������Բ���)�����������ɼнǵ�����ֵ������ֵ

real(8),save::strs_yld_0(30),strs_yld(800000)!����ͬ�Բ��ϵ�����ǿ��,ǰ��Ϊ�Բ��ʾ�������Ե�Ԫ��ʾ

real(8),save::strs_yld_110(30),strs_yld_220(30),strs_yld_330(30),strs_yld_120(30),strs_yld_130(30),strs_yld_230(30)  
real(8),save::strs_yld_11(800000),strs_yld_22(800000),strs_yld_33(800000),strs_yld_12(800000),strs_yld_13(800000),strs_yld_23(800000)
real(8),save::strs_yld_123(3,3,800000)
!�������Բ��ϵ�����ǿ��,ǰ��Ϊ�Բ��ʾ�������Ե�Ԫ��ʾ
real(8),save::f1_yield(800000)


real(8),save::strs_frct_0(30),strs_frct(800000)
real(8),save::strs_frct_123(3,3,800000)

real(8),save::strs_frct_110(30),strs_frct_220(30),strs_frct_330(30),strs_frct_120(30),strs_frct_130(30),strs_frct_230(30)
real(8),save::strs_frct_11(800000),strs_frct_22(800000),strs_frct_33(800000),strs_frct_12(800000),strs_frct_13(800000),strs_frct_23(800000)

integer(4),save::gas(800000)  !0�������ʲ������壬1��������Ϊ����

real(8),save::c0(30),s(30),gama(30)  !��̬���̲���
real(8),save::c0_iso(30),s_iso(30),gama_iso(30)
real(8),save::c0_aniso(30),s_aniso(30),gama_aniso(30)

real(8),save::h0(30)      !PUFF��̬���̲���
real(8),save::h0_iso(30) 
real(8),save::h0_aniso(30) 

real(8),save::cn(30)

real(8),save::es(30)         !������
real(8),save::es_iso(30)
real(8),save::es_aniso(30)

real(8),save::q1,q2    !ճ��ϵ��
real(8),save::qhg   !ɳ©ճ��ϵ�� ȡΪ0.1

real(8),save::vf

real(8),save::ax(30),ay(30),az(30),nx(30),ny(30),nz(30),axy(30),nxy(30),axz(30),nxz(30),ayz(30),nyz(30)
real(8),save::ax_aniso(30),ay_aniso(30),az_aniso(30),nx_aniso(30),ny_aniso(30),nz_aniso(30),axy_aniso(30),nxy_aniso(30),axz_aniso(30),nxz_aniso(30),ayz_aniso(30),nyz_aniso(30)
real(8),save::ax_iso(30),nx_iso(30)

real(8),save::rate(800000),strain_rate_eff(800000)
real(8),save::beta0(30),strain_rate0(30)        !Ӧ������Ч���ӵĲ���0.02183,��ʼ�ο�Ӧ���� ȡֵΪ0.001/s,Ӧ��������
real(8),save::beta0_iso(30),strain_rate0_iso(30)
real(8),save::beta0_aniso(30),strain_rate0_aniso(30)

integer(4),save::plastic_index(800000)   !���Ա���ָʾ���ӣ�Ϊ0ʱδ�������Ա��Σ�Ϊ1ʱ�������Ա���


real(8),save::u(3,800000),du(3,800000),v(3,800000),a(3,800000),v_mid(3,800000)  !��3nj�����ɶ��ϵĽڵ�λ��,�ڵ��ٶ�,�ڵ���ٶ�,�м�ʱ�̽ڵ��ٶ�
real(8),save::mass(800000)     !��nj�����ɶ��ϵ��������󣨶Խǻ���
real(8),save::x0(800000),y0(800000),z0(800000)   !��ʼʱÿ��Ԫ���ĵ���������
real(8),save::x(800000),y(800000),z(800000),v_cen_x(800000),v_cen_y(800000),v_cen_z(800000),l,h,m      !ÿ��Ԫ���ĵ��������꣬�ٶ�,��Ȼ����
!real(8),save::v_cen_mid_x(800000),v_cen_mid_y(800000)    !ÿ��Ԫ���Ĵ��������������м�ʱ�̵��ٶ�
!real(8),save::v_cen_mid_1(800000),v_cen_mid_2(800000)    !ÿ��Ԫ���Ĵ��������������м�ʱ�̵��ٶ�
! real(8),save::velocity_xyz(3,8,800000)
real(8),save::shape_n(8),shape_n0(8) !��״����,���ĵ����״����
real(8),save::dn_l(8),dn_h(8),dn_m(8),dn0_l(8),dn0_h(8),dn0_m(8)  !��״���������ĵ����״��������Ȼ����l,h,m�ĵ���
real(8),save::dn_lhm(8,3)
real(8),save::dn_xyz(8,3)
real(8),save::coord_xyz(3,8)
real(8),save::coord_xyz0(3,8)
real(8),save::dn_x(8),dn_y(8),dn_z(8)           !��״��������������x,y,z�ĵ���
real(8),save::jj(3,3),jj_inverse(3,3)       !(x,y,z)��(l,h,m)ƫ�������Ÿ��Ⱦ���,��������ʽ�������
real(8),save::jj0(3,3),jj0_inverse(3,3)    !(x0,y0,z0)��(l,h,m)ƫ�������Ÿ��Ⱦ���,��������ʽ�������
real(8),save::f11,f12,f13,f21,f22,f23,f31,f32,f33                  !(x,y,z)��(x0,y0,z0)��ƫ�������Ÿ��Ⱦ��󣬱����ݶȾ���
real(8),save::deform(3,3)
real(8),save::rotate_R(3,3)
real(8),save::rotate_original(3,3)


real(8),save::rho(800000),rho_old(800000),volume0(800000),volume(800000),volume_old(800000),dvolume(800000),volume_average(800000)
real(8),save::vol0(800000),vol(800000),vol_old(800000),dvol(800000),vol_average(800000)  !��Ԫ�ܶȣ���ʱ�������һʱ�̵����
real(8),save::e(800000),de(800000),p(800000),dp(800000),ee         !����Ԫ������,ƽ��Ӧ��,������

real(8),save::np(800000)

real(8),save::p_dst,p_dsp         !ƫӦ�乱��ѹ��������Ӧ�乱��ѹ��

real(8),save::f_ext(3,800000),f_int(3,800000),f_ext_old(3,800000),f_int_old(3,800000) 


real(8),save::stress_xx(800000),stress_yy(800000),stress_zz(800000),stress_xy(800000),stress_xz(800000),stress_yz(800000) !ÿ����Ԫ��l=0,s=0,m=0������Ӧ��
real(8),save::stress_xyz(3,3,80000)
real(8),save::stress_11(800000),stress_22(800000),stress_33(800000),stress_12(800000),stress_13(800000),stress_23(800000) !�����������У�ÿ����Ԫ��l=0,s=0,m=0������Ӧ��
real(8),save::stress_123(3,3,80000)
real(8),save::stress_d_11(800000),stress_d_22(800000),stress_d_33(800000),stress_d_12(800000),stress_d_13(800000),stress_d_23(800000)   !ÿ����Ԫ��l=0,s=0,m=0������ǰʱ�̵�ƫӦ��
real(8),save::stress_d_123(3,3,800000)

real(8),save::stress_node_xx(800000),stress_node_yy(800000),stress_node_zz(800000),stress_node_xy(800000),stress_node_xz(800000),stress_node_yz(800000),p_node(800000),e_node(800000) 
real(8),save::stress_node_xyz(3,3,80000)

real(8),save::diag_vector(3)
real(8),save::diag_matrix(3,3)
real(8),save::diag_sum
integer,save::number1,number2,number3
!��ƫ�ֽ�ר�ñ���

real(8),save::q(800000),q_old(800000),q_av(800000)   !�˹�ճ�����ǰʱ�̺���һʱ�̵�ֵ��ƽ��ֵ��
!==========================================================================================��ά��ɳ©����
real(8),save::hggama1(800000),hggama2(800000),hggama3(800000),hggama4(800000)  !ÿ����Ԫɳ©��ʸ���õ��ĸ���������ά������ĸ��Ƿ񹻣�����
real(8),save::f_hg(3,800000)     !ÿ���ڵ㴦��x,y�����ɳ©��

real(8),save::strain_xx(800000),strain_yy(800000),strain_zz(800000),strain_xy(800000),strain_xz(800000),strain_yz(800000)      !ÿ����Ԫ��l=0,s=0,m=0������Ӧ��,��Ч����Ӧ����
real(8),save::strain_xyz(3,3,80000)

real(8),save::strain_xx_rate(800000),strain_yy_rate(800000),strain_zz_rate(800000),strain_xy_rate(800000),strain_xz_rate(800000),strain_yz_rate(800000)   !Ӧ���ʷ���
real(8),save::strain_xyz_rate(3,3,800000)
real(8),save::strain_11(800000),strain_22(800000),strain_33(800000),strain_12(800000),strain_13(800000),strain_23(800000)      !�����������У�ÿ����Ԫ��l=0,s=0,m=0������Ӧ��
real(8),save::strain_123(3,3,800000)

real(8),save::strain_11_rate(800000),strain_22_rate(800000),strain_33_rate(800000),strain_12_rate(800000),strain_13_rate(800000),strain_23_rate(800000)   !��������ϵ��Ӧ���ʷ���
!real(8),save::strain_123_rate(3,3,800000)
real(8),save::strain_p_11(800000),strain_p_22(800000),strain_p_33(800000),strain_p_12(800000),strain_p_13(800000),strain_p_23(800000)      !�����������У�ÿ����Ԫ��l=0,s=0,m=0����������Ӧ��
real(8),save::strain_p_123(3,3,800000)

real(8),save::strain_p_eff(800000),strain_p_rate_eff(800000)  !��Ч����Ӧ��,��Ч����Ӧ����������Ч����Ӧ����

!real(8),save::peak_axis(800000)  !�����ϸ���Ԫ��Ӧ����ֵ

real(8),save::rad_energy_total  !����ͨ���ܶ�
integer(4),save::rad_time_mode  !ʱ����ģʽ��0Ϊ�����ף�1Ϊ��������
real(8),save::rad_time_width    !ʱ���׵Ŀ�ȣ���������Ϊ��߿�
integer(4),save::rad_energy_num !ʱ�����=ʱ���׿��/dt
real(8),save::rad_energy(8000000)   !�����ܶ�(energy_numά)

integer(4),parameter::ken=5000
real(8),save::wl_short,wl_long
real(8),save::wl(5000)
real(8),save::wer(5000)

integer(4),save::n_blackbody  !�����׵ĸ��������ȡ20��
real(8),save::t_kev(20)       !��������ֵĺ����¶ȣ���λΪkeV
real(8),save::t_k(20)       !��������ֵĺ����¶ȣ���λΪk(��)
real(8),save::b_energy(20)    !�����������ռ�������������ݶ�
real(8),save::f_blackbody_each(20)  !��������ֵ����ʿ˺���
real(8),save::a_planck(20)    !�������׷���ǿ�ȵı���ϵ��

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
    integer(4) element_num    !�����Ԫ��������<=10

	integer(4) element_z(max_element_num)    !Ԫ�ص�ԭ��������Ԫ�����಻����10��
	real(8) element_wp(max_element_num)     !��Ԫ�ص������ٷֱ�

	real(8) mu(ken)
end type layer_xray

type(layer_xray), save :: xlayer(max_layer)

real(8),save::w_int,w_ext,w_kin,w_orgn,w_err     !���ܣ����ܣ����ܣ���ʼ�������������ٷֱ�
!���������ر���
real,save::core(3,1)
real,save::volume_s(6)






end module vars   