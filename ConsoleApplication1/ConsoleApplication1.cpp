// ConsoleApplication1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <math.h>
#include <limits>
#include <vector>
using namespace std;

template<class Type>
void print_matrix(Type** matrix, int size_x, int size_y) {
	for (int i = 0; i < size_x; i++) {
		for (int j = 0; j < size_y; j++) {
			cout << matrix[i][j];
			if (j != size_y - 1) {
				cout << "\t";
			}
		}
		cout << "\n";
	}
}

float distance(float* p1, float* p2) {
	return sqrt(pow(p2[0] - p1[0], 2) + pow(p2[1] - p1[1], 2));
}

//генерация матрицы смежности графа на основе расстояний между n случайно сгенерированнными  в единичном квадрате  точками
float** generateFlatGraph(int n) {
	//параметры генератора случайных чисел
	const float accuracy = 0.01;
	const int divider = (int) 1 / accuracy;

	//генерация координат точек
	float** point_coord = new float*[n];
	for (int i = 0; i < n; i++) {
		point_coord[i] = new float[2];
		point_coord[i][0] = (float)(rand() % (divider + 1)) / divider;
		point_coord[i][1] = (float)(rand() % (divider + 1)) / divider;
	}

	//вычисление матрицы смежности графа на основе расстояний между точками
	float** adj_matrix = new float*[n];
	for (int i = 0; i < n; i++) {
		adj_matrix[i] = new float[n];
		for (int j = 0; j < n; j++) {
			if (i != j) {
				adj_matrix[i][j] = distance(point_coord[i], point_coord[j]);
			}
			else
			{
				adj_matrix[i][j] = numeric_limits<float>::max();
			}
		}
	}

	return adj_matrix;
}

template <class Type>
float min(Type* mas, int size) {
	Type min_val = numeric_limits<Type>::max();
	for (int j = 0; j < size; j++) {
		if (mas[j] < min_val) {
			min_val = mas[j];
		}
	}
	return min_val;
}

template <class Type>
float reduce(Type* mas, int size, Type val) {
	for (int i = 0; i < size; i++) {
		mas[i] = mas[i] - val;
	}
	return val;
}

template <class Type>
float min_column(Type** matrix, int column, int size) {
	Type min_val = numeric_limits<Type>::max();
	for (int j = 0; j < size; j++) {
		if (matrix[j][column] < min_val) {
			min_val = matrix[j][column];
		}
	}
	return min_val;
}

template <class Type>
float reduce_column(Type** matrix, int column, int size, Type val) {
	for (int i = 0; i < size; i++) {
		matrix[i][column] = matrix[i][column] - val;
	}
	return val;
}

template <class Type>
float reduce_matrix(Type** matrix, int size_x, int size_y) {
	Type val = 0;
	for (int i = 0; i < size_x; i++) {
		val += reduce(matrix[i], size_x, min(matrix[i], size_x));
	}
	for (int i = 0; i < size_y; i++) {
		val += reduce_column(matrix, i, size_y, min_column(matrix, i, size_y));
	}
	return val;
}

template <class Type>
vector<int*> matrix_null_element_positions(Type** matrix, int size) {
	vector<int*> null_element_positions;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (matrix[i][j] == 0) {
				int* pos = new int[2]{i,j};
				null_element_positions.push_back(pos);
			}
		}
	}
	return null_element_positions;
}

int main()
{
	const int n = 5;
	float** matrix = generateFlatGraph(n);
	/*
	reduce_matrix(matrix, n, n);
	print_matrix(matrix, n, n);
	vector<int*> null_pos = matrix_null_element_positions(matrix, n);
	while (!null_pos.empty()) {
		cout << null_pos.back()[0] << "\t" << null_pos.back()[1] << "\n";
		null_pos.pop_back();
	}
	*/
}
