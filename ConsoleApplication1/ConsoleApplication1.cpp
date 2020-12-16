// ConsoleApplication1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <math.h>
#include <limits>
using namespace std;

//вывод матрицы
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

//расстояние между точками p1 и p2
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

	print_matrix(point_coord, n, 2);
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


int main()
{
	print_matrix<float>(generateFlatGraph(5), 5, 5);
}
