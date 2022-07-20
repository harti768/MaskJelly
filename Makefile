#Compiler
CC = g++
FLAGS = -Wall -g -pthread

#Source folder
SRC = $(realpath $(CURDIR))/src/

all: maskjelly

maskjelly: mask_mers.o command_line_parser.o
	$(CC) $(FLAGS) $(INC) -o maskjelly mask_mers.o command_line_parser.o

mask_mers.o: $(SRC)mask_mers.cpp
	$(CC) $(FLAGS) $(INC) -c $(SRC)mask_mers.cpp

command_line_parser.o: $(SRC)command_line_parser.cpp $(SRC)command_line_parser.hpp
	$(CC) $(FLAGS) $(INC) -c $(SRC)command_line_parser.cpp

clean:
	rm -f mask_mers.o command_line_parser.o maskjelly
