TARGET = mix

CC = g++
CFLAGS = -O3 -std=c++11 -DDOUBLE

SRC = \
	./src/main.cpp

INC = ./src/

default:
	make clean
	make $(TARGET)

$(TARGET):
	$(CC) $(CFLAGS) $(SRC) -I$(INC) -o $(TARGET)

clean:
	rm -f $(TARGET) 
