program bin
	use SHTOOLS
	implicit none
	character ::		filein*80,b(12000)*8,filetopo*80,table(100)*80,filefa*80,filedegree*80,fileerror*80
	real*8 ::		d(100),e(180,360),f(6480,10),g(180,360),mpr,ss,&
				pi
	real*8,allocatable ::	coef_c(:,:),coef_s(:,:),m_covar(:,:),a_coef(:,:),a_covar(:),m_covar_seq(:,:),a_coef_seq(:,:),&
				grid(:,:),grid_out(:,:),pl(:,:),d_cnm_l(:),d_snm_l(:),d_g_l(:),d_gp_l(:),grid_error(:,:),grid_degree(:,:),&
				topoc(:,:,:),grid_topo(:,:)
	integer*8 ::		i,j,img(898,8),ff(1025,1025),ff2(1024),ii(100),nmax,n_coef,n_var,n_degree,p,n_cnm,&
				i_cnm,i_snm,l,l_cnm,m,n,ll,selectdegree,nline_end,n_error_degree,l_error
	integer ::		lmaxt 
	real*8 			ae, gm, gmsig, reflon, reflat,grid_int,theta,lat,lon,fa,fa_l,d_cnm,d_snm,error_g,kaula_constant,&
				kaula_k,timein,timeout,para(8),unuse
	integer*4 		lmax,mmax,inorm, nvar, nline, jend
	pi = acos(-1d0)
	nmax = 1d5
	

	!read(10,rec=1)ae, gm, gmsig, lmax, mmax, inorm, nvar, reflon, reflat
	selectdegree = 100
	select case(selectdegree)
		case(100)
		filein = "GGMES_100V08_SHB.DAT"
		kaula_constant = 3d-5
		case(50)
		filein = "ggmes_50v06_shb.dat"
		kaula_constant = 1.25d-5
		case(20)
		filein = "ggmes_20v02_shb.dat"
		kaula_constant = 1.25d-5
		case(95)
		filein = "jgmro_095a_shb.dat"
		kaula_constant = 1.25d-5
	end select
	filetopo = "gtmes_150v05_sha.tab"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! reading header files
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (selectdegree==95) then
	open (12, file=filein, status ='old', form = 'unformatted',access='DIRECT',convert='little_endian', RECL=512)
	else
	open (12, file=filein, status ='old', form = 'unformatted',access='DIRECT',convert='big_endian', RECL=512)
	end if
	read(12,rec=1) ae, gm, gmsig, lmax, mmax, inorm, nvar, reflon, reflat

	print*, "average radius:",ae,"GM:", gm,"GM error:", gmsig
	print*,"normal gravity",gm/(ae**2)*1d5*1d3
	print*, "maxium degree",lmax, "maxium order", mmax, "normaliz state", inorm, "number of parameters",nvar

	print*, "reference longitute",reflon, "reference latitute",reflat
	nline = (nvar/64) + 1
	print*,"number of line",nline
	jend = nvar - (nvar/64)*64
	print*,jend
	n_degree = lmax
	!n_degree = 20
	j=0
	do i =2,n_degree
		j=j+i+1
	end do
	print*,"Cnm number",j
	n_cnm=j
	do i =2,n_degree
		j=j+i
	end do
	print*,"coeffient number",j
	n_coef = j
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate covariance parameters number
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	j=0
	do i =1,nvar
		j=j+i
	end do
	print*,"covariance number",j
	n_var = j
	!print*,n_var/64 + 1
	print*,"reading header over"

	print*,"maxium degree",int(n_degree,4),"number of Cnm",int(n_cnm,4)
	print*,"number of Cnm+Snm",int(n_coef,4),"number of parameters",int(nvar,4)
	print*,"therotical variance number",int(n_var,4),"nvar**2",int(nvar*nvar,4)
	print*,"reading variance from line",nline *2 + 2,"to line ",nline *2 + 2 + int(n_var/64)
	allocate(coef_c(120, 120))
	allocate(coef_s(120, 120))
	allocate(m_covar(12000, 12000))
	allocate(m_covar_seq(12000, 12000))
	allocate(grid(12000, 12000))
	allocate(grid_error(12000, 12000))
	allocate(grid_degree(12000, 12000))
	allocate(grid_topo(12000, 12000))
	allocate(grid_out((n_degree+1)*(n_degree*2+1), 3))
	allocate(a_coef(12000,4))
	allocate(a_coef_seq(12000,4))
	allocate(d_cnm_l(12000))
	allocate(d_snm_l(12000))
	allocate(d_g_l(12000))
	allocate(d_gp_l(12000))
	allocate(a_covar(60000000))
	allocate(pl(600000,4))
	allocate(topoc(2,720,720))
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! reading topography file
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call SHRead(filetopo, topoc, lmaxt ,header=para(1:2))
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! reading coeffient name table
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i = 2,nline+1
	read(12,rec=i) b((i-2)*64+1:(i-2)*64+64)
	end do

	do i = 1,nline*64

	end do

	close(12)
	!go to 97
	if (selectdegree==95) then
	open (12, file=filein, status ='old', form = 'unformatted',access='DIRECT',convert='little_endian', RECL=512)
	else
	open (12, file=filein, status ='old', form = 'unformatted',access='DIRECT',convert='big_endian', RECL=512)
	end if
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! reading coeffient value table
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i = nline + 2,nline *2 + 1
	read(12,rec=i) a_coef((i-2-nline)*64+1:(i-2-nline)*64+64,3)
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! test parameters order
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!do i = 1,n_coef	+ 10
	do i = nvar-50,nvar+50
	!print*,b(i),a_coef(i,3),i
	end do


	close(12)
	print*,"reading coeffient table over"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! assign degree and order to the coeffient : Cnm and Snm
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i = 0,n_degree
		do j = 0,i
			if ((i)*(i+1)/2 + j -2 >0) then
		!print*,i+1,j+1,(i)*(i+1)/2 + j -2
		coef_c(i+1,j+1) = a_coef((i)*(i+1)/2 + j -2,3)
		a_coef((i)*(i+1)/2 + j -2,1) = i
		a_coef((i)*(i+1)/2 + j -2,2) = j
		end if
		!print*,i,j, coef_c(i+1,j+1)
		end do
	end do
	do i = 0,n_degree
		do j = 0,i
		!print*,i+1,j+1,(i)*(i+1)/2 + j -3
		if (j ==0) then
		!print*,i+1,j+1,0
		coef_s(i+1,j+1) = 0d0
		else if ((i)*(i+1)/2 + j -1 -i >0) then
		!print*,i+1,j+1,(i)*(i+1)/2 + j -1 -i
		a_coef((i)*(i+1)/2 + j - 2 -i + n_cnm+1,1)=i
		a_coef((i)*(i+1)/2 + j - 2 -i + n_cnm+1,2)=j
		coef_s(i+1,j+1) = a_coef((i)*(i+1)/2 + j -2 -i + n_cnm+1,3)
		end if
		!print*,i,j, coef_c(i+1,j+1)

		end do
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! test the degree and order
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i = 0,n_degree
		do j = 0,i
		!print*,i,j, coef_c(i+1,j+1), coef_s(i+1,j+1)

		end do
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! test the sequence
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i = 1,nline*64
		!print*,b(i),a_coef(i,1),a_coef(i,2),a_coef(i,3),i
	end do


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! reading covariance table
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!go to 97
	if (selectdegree==95) then
	open (12, file=filein, status ='old', form = 'unformatted',access='DIRECT',convert='little_endian', RECL=512)
	else
	open (12, file=filein, status ='old', form = 'unformatted',access='DIRECT',convert='big_endian', RECL=512)
	end if
	select case(selectdegree)
	case(20)
		nline_end = nline *2 + 2 + int(n_var/64)
	case(50)
		nline_end = nline *2 + 2 + int(n_var/64)
	case(100)
		nline_end = nline *2 + 2 + int(n_var/64)
	case(95)
		nline_end = nline *2 + 2 + int(n_var/64)
	end select
	print*,"nline_end",int(nline *2 + 2,4),int(nline_end,4)
	do i = nline *2 + 2,nline_end
		read(12,rec = i) a_covar((i - nline*2 - 2)*64 + 1:(i-nline*2- 2)*64 + 64)

	end do
	print*,"reading variance from line",nline *2 + 2,"to line ",int(nline_end,4)
	close(12)

	do i = n_var - 5,n_var + 3
	!print*,a_covar(i),i
	end do
	do i = 1,80
	!print*,a_covar(i),i,a_covar(i+nvar-1),i	+ nvar - 1
	end do
	do i = nvar-3,nvar+5
	!print*,a_covar(i),i
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! transit the covariance to matrix
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do i = 1,nvar
		if (selectdegree ==95) then
		m_covar(1:i,i) = a_covar((i-1)*i/2+1:(i+2)*(i-1)/2+1)
		m_covar(i,1:i) = a_covar((i-1)*i/2+1:(i+2)*(i-1)/2+1)
		else
		m_covar(i,i:nvar) = a_covar((i-1)*(nvar + 1 ) + 1-i*(i-1)/2:i*nvar-i*(i-1)/2)
		m_covar(i:nvar,i) = a_covar((i-1)*(nvar + 1 ) + 1-i*(i-1)/2:i*nvar-i*(i-1)/2)
		!print*,(i-1)*438+1-i*(i-1)/2,i*437-i*(i-1)/2,i
		end if
	end do
	go to 90
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! test the matrix order
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	print*,"full covariance"
	write(*,*)
	print*,"first 5 "

		do j = 1,5
		print*,real(m_covar(j,1:5),4)
		write(*,*)
		end do

	write(*,*)
	print*,"lase 5 "

		do j = nvar - 4,nvar
		print*,real(m_covar(j,nvar - 4:nvar),4)
		write(*,*)
		end do

	write(*,*)
	90 continue
	do i = 1, 100
		!print*,b(i),real(a_coef(i,3),4),m_covar(i,i),int(i,4)
	end do
	!print*,m_covar(1,1)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! generate index of each degree and order
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!go to 97
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!
	! the kaula law
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!
	n_error_degree = n_degree
	kaula_k = gm/ae**2*kaula_constant*sqrt(2d0)*1d3*1d5
	!print*,kaula_k
	do i = 1,n_error_degree
		!print*,kaula_f(i,kaula_constant,gm,ae)*1d8,kaula_k/sqrt(dble(i))
	end do
	!n_error_degree = 100
	grid_degree(:,:)=0d0
	
	call cpu_time(timein)
	n_degree = 100
	!go to 89
	do l_error = 1,n_error_degree

		i_cnm = 0
		i_snm = 0
		m_covar_seq(:,:) = 0d0
		a_coef_seq(:,:) = 0d0
		do i = 2,l_error
			i_cnm = i_cnm + i + 1
			i_snm = i_snm + i
		end do
		l_cnm = i_cnm+i_snm
		print*,int(l_error,4),int(i_cnm,4),int(n_cnm,4),int(n_cnm+i_snm,4),int(i_snm,4),int(i_cnm+i_snm,4)
		m_covar_seq(1:i_cnm,1:i_cnm) = m_covar(1:i_cnm,1:i_cnm)
		m_covar_seq(1:i_cnm,i_cnm + 1:i_cnm+i_snm) = m_covar(1:i_cnm,n_cnm + 1:n_cnm+i_snm)
		m_covar_seq(i_cnm + 1:i_cnm+i_snm,1:i_cnm) = m_covar(n_cnm + 1:n_cnm+i_snm,1:i_cnm)
		m_covar_seq(i_cnm + 1:i_cnm+i_snm,i_cnm + 1:i_cnm+i_snm) = m_covar(n_cnm + 1:n_cnm+i_snm,n_cnm + 1:n_cnm+i_snm)
		a_coef_seq(1:i_cnm,1:3) = a_coef(1:i_cnm,1:3)
		a_coef_seq(i_cnm + 1:i_cnm+i_snm,1:3) = a_coef(n_cnm + 1:n_cnm+i_snm,1:3)
		
	grid_int = 180d0/dble(n_degree)
	do i = n_degree + 1,1,-1
	!do i = 1, n_degree + 1
		do j = n_degree*2 + 1,1,-1
		!do j = 1, n_degree*2 + 1
			if (grid_degree(i,j) == 0) then
			lat = (dble(i-1)*grid_int)*pi/180d0
			lon = (180d0 -dble(j-1)*grid_int)*pi/180d0
			fa = 0d0
			d_cnm = 0d0
			d_snm = 0d0
			error_g = 0d0
			call PlmBar(pl(:,3), int(200,4), cos(lat), csphase = 1, cnorm = 0)
			!call PLegendreA(pl(:,3), int(n_degree,4), cos(lon), csphase = 1)
			do l = 1, i_cnm
				fa_l = a_coef_seq(l,3)*cos(a_coef_seq(l,2)*lon)*&
				pl(index_pl(int(a_coef_seq(l,1),8),int(a_coef_seq(l,2),8)),3)*(a_coef_seq(l,1) + 1d0)
				fa = fa + fa_l
				d_g_l(l) = cos(a_coef_seq(l,2)*lon)*&
				pl(index_pl(int(a_coef_seq(l,1),8),int(a_coef_seq(l,2),8)),3)*(a_coef_seq(l,1) + 1d0)*&
				gm/(ae**2)*1d5*1d3

			end do
			do l = i_cnm + 1, i_cnm + i_snm
				fa_l = a_coef_seq(l,3)*sin(a_coef_seq(l,2)*lon)*&
				pl(index_pl(int(a_coef_seq(l,1),8),int(a_coef_seq(l,2),8)),3)*(a_coef_seq(l,1) + 1d0)
				fa = fa + fa_l
				d_g_l(l) = sin(a_coef_seq(l,2)*lon)*&
				pl(index_pl(int(a_coef_seq(l,1),8),int(a_coef_seq(l,2),8)),3)*(a_coef_seq(l,1) + 1d0)*&
				gm/(ae**2)*1d5*1d3
			end do
			d_gp_l(:) = 0d0
			do l = 1, i_cnm + i_snm
				do ll = 1, i_cnm + i_snm
				d_gp_l(l) = d_gp_l(l) + d_g_l(ll)*m_covar_seq(l,ll)
				end do
				d_gp_l(l) = d_gp_l(l)*d_g_l(l)
			end do
			error_g = sqrt(sum(d_gp_l(1:i_cnm + i_snm)))
			if(error_g > kaula_k/sqrt(dble(l_error)) .or. l_error == n_error_degree) then
			fa = fa*gm/(ae**2)*1d5*1d3
			grid_error(i,j) = error_g
			grid(i,j) = fa
			grid_degree(i,j) = l_error
			end if
			end if
			!print*,"fa",real(dble(i-1)*grid_int,4),real(180d0 -dble(j-1)*grid_int,4),fa,error_g
		end do
	end do
	end do
	89 continue
	call cpu_time(timeout)
	print*, "time (sec) = ", timeout-timein
	do i = 1,n_degree + 1
	!do j = 1,n_degree*2 + 1
		!print*,grid(i,1)
		!print*,grid(1,j)
	end do
	open(12, file = "grid_fa_100_dsm.dat")
	open(13, file = "grid_error_100_dsm.dat")
	open(14, file = "grid_degree_100_dsm.dat")
	do i = 1, n_degree + 1
		do j = 1, n_degree*2 + 1
			!read(12,*) unuse,unuse, grid(i,j)
			!read(13,*) unuse,unuse, grid_error(i,j)
			!read(14,*) unuse,unuse, grid_degree(i,j)
		end do
	end do
	close(12)
	close(13)
	close(14)
	open(12, file = "grid_fa.dat")
	open(13, file = "grid_error.dat")
	open(14, file = "grid_degree.dat")
	do i = 1, n_degree + 1
		do j = 1, n_degree*2 + 1
			write(12,*) real(180d0 - dble(j-1)*grid_int,4),real(90d0 - dble(i-1)*grid_int, 4), grid(i,j)
			write(13,*) real(180d0 - dble(j-1)*grid_int,4),real(90d0 - dble(i-1)*grid_int, 4), grid_error(i,j)
			write(14,*) real(180d0 - dble(j-1)*grid_int,4),real(90d0 - dble(i-1)*grid_int, 4), grid_degree(i,j)
		end do
	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate coeffient number
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	j=0
	do i =2,20
		j=j+i+1
	end do
	print*,"Cnm number",j

	do i =2,20
		j=j+i
	end do
	print*,"coeffient number",j

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate covariance parameters number
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	j=0
	do i =1,437
		j=j+i
	end do
	print*,"covariance number",j



	print*,4.869d24*6.67d-11/6050d3/6050d3*1.2d-5*sqrt(2d0)*1d5
	theta = 30d0*pi/180d0
	call PLegendreA(pl(:,3), int(n_degree,4), cos(theta))
	i=1
	do n = 0,n_degree
		do m = 0,n
			pl(i,1) = n
			pl(i,2) = m
			i = i + 1
		end do
	end do
	go to 99
	print*,i
	print*,1
	print*,cos(theta)
	print*,sin(theta)
	print*,75d-2*cos(theta*2) + 25d-2
	print*,3d0*sin(theta)*cos(theta)
	print*,3*sin(theta)**2
	print*,5d0/8d0*cos(theta*3) + 3d0/8d0*cos(theta)
	print*,sin(theta)*(15d0/2d0*cos(theta)**2 - 3d0/2d0)
	print*,15*sin(theta)**2*cos(theta)
	print*,15*sin(theta)**3
	write(*,*)
	do i = 1,10
	print*,index_pl(int(pl(i,1),8),int(pl(i,2),8)),pl(i,1),pl(i,2),pl(i,3)
	end do
	99 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














	100 format(D17.8,1X,D17.8,1X,F23.8)
	101 format(I8,1X,I8,1X,I8,1X,I8)
	102 format(Z16,Z16,Z16,Z8,Z8,Z8,Z8,Z16,Z16)
	103 format(Z16,1X,Z16,1X,Z16,1X,Z8,1X,Z8,1X,Z8,1X,Z8,1X,Z16,1X,Z16)
	104 format(B64,B64,B64)
	105 format(Z16,1X,Z16,1X,Z16)
contains
	function index_pl(n,m)
		integer*8 m,n,i,j,index_pl
		index_pl = 1
		do i = 0,n-1
		do j = 0,i

			index_pl = index_pl + 1
		end do
		end do
		index_pl = index_pl + m
	end function index_pl
	function kaula_f(n,k,gm,r)
		integer*8 n
		real*8 kaula_f,k,gm,r
		kaula_f = gm/r**2*k*sqrt(2d0/dble(n))
	end function









end program bin
