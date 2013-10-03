all: cole_app

cole_app: multi_elliptical_orbit_tracking_complicated_fit.o
	c++ multi_elliptical_orbit_tracking_complicated_fit.o -o cole_app 

multi_elliptical_orbit_tracking_complicated_fit.o: multi_elliptical_orbit_tracking_complicated_fit.cpp nr.h nrutil.h
	c++ -c multi_elliptical_orbit_tracking_complicated_fit.cpp

clean:
	rm -rf *o
