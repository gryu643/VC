program ppl_bulk
	implicit none

	integer,parameter :: row=3,col=3,B_col=1
	integer i,j

	complex A(row,col)
	complex B(row,B_col)
	complex C(row,B_col)

	do i=1, col
		do j=1, col
			if(i==2) then
				A(i,j)=(2.0, 3.0)
			else
				A(i,j)=(1.0, 1.0)
			end if
		end do
	end do

	do i=1, row
		B(i, B_col) = (1.0, 1.0)
		C(i, B_col) = (0.0, 0.0)
	end do

	call Multiply(A,B,C,row,col,row,B_col)

	do i=1, 3
		print *, C(i,B_col)
	end do

end program 

subroutine Multiply(A,B,C,A_row,A_col,B_row,B_col)
	integer A_row,A_col,B_row,B_col
	complex A(A_row,A_col), B(B_row,B_col), C(A_row,B_col)

	if(A_col.ne.B_row) then
		print *, "can't caluculate (Multiply)"
		stop
	end if

	do i=1, A_row
		do j=1, B_col
			do k=1, A_col
				C(i,j) = C(i,j)+A(i,k)*B(k,j)
			end do
		end do
	end do

end subroutine