#include <iostream>
#include <math.h>
#include <limits>
#include <vector>
#include <type_traits>
#include <fstream>
#include <string>
using namespace std;

template<class Type>
void print_matrix(Type** matrix, int size_x, int size_y) {
	for (int i = 0; i < size_x; i++) {
		for (int j = 0; j < size_y; j++) {
			cout << matrix[i][j];
			if (j != size_y - 1) {
				cout << "\t\t";
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
Type min(Type* mas, int size) {
	Type min_val = numeric_limits<Type>::max();
	for (int j = 0; j < size; j++) {
		if (mas[j] < min_val) {
			min_val = mas[j];
		}
	}
	return min_val;
}

template <class Type>
Type reduce(Type* mas, int size, Type val) {
	for (int i = 0; i < size; i++) {
		mas[i] = mas[i] - val;
	}
	return val;
}

template <class Type>
Type min_column(Type** matrix, int column, int size) {
	Type min_val = numeric_limits<Type>::max();
	for (int j = 0; j < size; j++) {
		if (matrix[j][column] < min_val) {
			min_val = matrix[j][column];
		}
	}
	return min_val;
}

template <class Type>
Type reduce_column(Type** matrix, int column, int size, Type val) {
	for (int i = 0; i < size; i++) {
		matrix[i][column] = matrix[i][column] - val;
	}
	return val;
}

template <class Type>
Type reduce_matrix(Type** matrix, int size_x, int size_y) {
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

template <class Type>
Type zero_degree(Type** matrix, int* pos, int n) {
	Type t = matrix[pos[0]][pos[1]];
	matrix[pos[0]][pos[1]] = numeric_limits<Type>::max();
	Type min_val = min(matrix[pos[0]], n) + min_column(matrix, pos[1], n);
	matrix[pos[0]][pos[1]] = t;
	return min_val;
}

template <class Type>
Type max_zero_degree_null(Type** matrix, int* pos, int n) {
	Type max = numeric_limits<Type>::min();
	vector<int*> null_positions = matrix_null_element_positions(matrix, n);
	while (!null_positions.empty()) {
		Type val = zero_degree(matrix, null_positions.back(), n);
		if (val > max) {
			max = val;
			pos[0] = null_positions.back()[0];
			pos[1] = null_positions.back()[1];
			null_positions.back();
		}
		null_positions.pop_back();
	}
	return max;
}

template <class Type>
Type** matrix_from_file(string path) {
	ifstream in(path);
	string line;
	vector<vector<Type>> matrix_lines;
	if (in.is_open())
	{
		while (getline(in, line))
		{
			vector<Type> matrix_line;
			int last_space_pos = 0;
			int space_pos = line.find('\t', last_space_pos);
			do
			{
				string elem = line.substr(last_space_pos, space_pos != -1 ? space_pos - last_space_pos : line.length());
				matrix_line.push_back(elem != "-" ?
					is_same<Type, int>::value ? stoi(elem) :
					is_same<Type, float>::value ? stof(elem) :
					numeric_limits<Type>::max() - 100 :
					numeric_limits<Type>::max() - 100
				);
				last_space_pos = space_pos + 1;
				space_pos = line.find('\t', last_space_pos);
			} while (last_space_pos != 0);
			matrix_lines.push_back(matrix_line);
		}
		in.close();
		Type** matrix = new Type*[matrix_lines.size()];
		for (int i = 0; i < matrix_lines.size(); i++) {
			matrix[i] = new Type[matrix_lines.at(i).size()];
			for (int j = 0; j < matrix_lines.at(i).size(); j++) {
				matrix[i][j] = matrix_lines.at(i).at(j);
			}
		}
		return matrix;
	}
	return NULL;
}

template <class Type>
struct tree_node {
	int* edge;
	Type val;
};

int main()
{
	const int n = 5;
	//float** matrix = generateFlatGraph(n);

	int** matrix = matrix_from_file<int>("input.txt");
	
}
