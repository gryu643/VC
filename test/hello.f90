program hello
	integer,parameter :: x=4,y=12
	integer z

	call sum(x,y,z)

	print *, z

end program

subroutine sum(a,b,c)
	integer a,b,c
	c = a + b
end subroutine