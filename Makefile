TARGET = mix

default:
	make clean
	make cpp

cpp:
	make --directory=src_cpp
	mv src_cpp/mix $(TARGET)

fortran:
	make --directory=src_fortran
	mv src_fortran/mixing_test $(TARGET)
	
clean:
	rm -f $(TARGET)
