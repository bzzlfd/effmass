all: convert_wg2rho.x 
convert_wg2rho.x:
	ifort convert_wg2rho.f90 cfft.f90 cfftd.f90 -o $@

clean:
	rm -f *.x
