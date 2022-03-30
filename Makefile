TARGET = mix

CC = g++
CFLAGS = -O3 -std=c++11 -DDOUBLE

SRC = \
	./src_c/main.cpp

INC = ./src_c/

default:
	make clean
	make $(TARGET)

$(TARGET):
	$(CC) $(CFLAGS) $(SRC) -I$(INC) -o $(TARGET)

clean:
	rm -f $(TARGET) 
