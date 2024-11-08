CNDO:
	cd src && g++ -std=c++14 -o fock Fock.cpp -I/opt/homebrew/include -L/opt/homebrew/lib -larmadillo && ./fock