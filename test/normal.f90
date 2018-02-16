program seikibunpu
	implicit none
	integer,parameter :: trial=100000
	double precision,parameter :: stride=0.1
	integer seedsize,i,j
	integer,allocatable :: seed1(:)
	integer,allocatable :: seed2(:)
	double precision :: out=0.0
	double precision :: result(100)=0.0

	!seedサイズの取得・動的配列の再宣言
	call random_seed(size=seedsize)
	allocate(seed1(seedsize))
	allocate(seed2(seedsize))

	!初期シードを設定（独立になるようseed配列を書き換える）
	do i=1, seedsize
		seed1(i) = 1
	end do
	do i=1, seedsize
		seed2(i) = 1000000
	end do

	do i=1, trial
		out = normal()

		do j=1, 100
			if((out.ge.(-5.0+stride*j)).and.(out.lt.(-5.0+stride*(j+1)))) then
				result(j) = result(j) + 1.0
			end if
		end do
	end do

	do i=1, 100
		write(*, fmt='(F6.2)', advance='no') -5.0+stride*i
		write(*, fmt='(a)', advance='no') ": "
		do j=1, int(result(i)/result(50)*50)
			write(*, fmt='(a)', advance='no') "*"
		end do
		print *, " "
	end do

!	call random()

contains
	subroutine random()
		double precision :: x=0.0
		double precision :: y=0.0
		integer seedsize,i
		integer,allocatable :: seed(:)

		call random_seed(size=seedsize)
		allocate(seed(seedsize))

		!シード１
		do i=1, 10
			call random_number(x)
			print *, i, x
			if(i.eq.5) then
				call random_seed(get=seed)
			end if
		end do

		!シード変更
		call random_seed(put=seed)

		do i=1, 5
			call random_number(x)
			print *, i, x
		end do
	end subroutine

	function normal()
		double precision :: m=0.0
		double precision :: var=1.0
		double precision :: x=0.0
		double precision :: y=0.0
		double precision :: pi=dacos(dble(-1))
		double precision :: r1=0.0
		double precision :: r2=0.0
		double precision normal

		call random_seed(put=seed1)
		call random_number(r1)
		call random_seed(get=seed1)

		call random_seed(put=seed2)
		call random_number(r2)
		call random_seed(get=seed2)


		x = sqrt(-2.0*dlog(r1))*dcos(2.0*pi*r2)
		y = sqrt(-2.0*dlog(r1))*dsin(2.0*pi*r2)

		normal = x
	end function

end program

