CXXFLAGS=-std=c++11 -Wpedantic -Wall -Wextra -Werror -Wzero-as-null-pointer-constant -Weffc++
CFLAGS=-std=c99 -pedantic -Wall -Wextra -Werror -ansi -Wwrite-strings -Weffc++

main: AVL.cpp
	echo     ------------compiling AVL.ccp to create executable program AVL----------------
	g++ -ggdb -std=c++11 AVL.cpp -o AVL

clean:
	echo ----------------removing executable program AVL--------------------
	/bin/rm AVL
