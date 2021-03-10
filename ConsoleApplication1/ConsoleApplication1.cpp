#include <iostream>
#include <math.h>
#include <limits>
#include <vector>
#include <type_traits>
#include <fstream>
#include <string>
using namespace std;

int transformIndex(int source_index, vector<int> index_list);

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

template<class Type> 
Type** clone_matrix(Type** matrix, int size_x, int size_y) {
	Type** new_matrix = new Type*[size_x];
	for (int i = 0; i < size_x; i++) {
		new_matrix[i] = new Type[size_y];
		for (int j = 0; j < size_y; j++) {
			new_matrix[i][j] = matrix[i][j];
		}
	}
	return new_matrix;
}

template <class Type>
struct numbered_matrix {
	Type** matrix;
	vector<int> out_indexes;
	vector<int> in_indexes;
};

template <class Type>
numbered_matrix<Type> remove_edge(numbered_matrix<Type> n_matrix, int* edge) {
	vector<int> new_out_indexes(n_matrix.out_indexes);
	vector<int> new_in_indexes(n_matrix.in_indexes);
	new_out_indexes.erase(new_out_indexes.begin() + transformIndex(edge[0], new_out_indexes));
	new_in_indexes.erase(new_in_indexes.begin() + transformIndex(edge[1], new_in_indexes));
	Type** new_matrix = new Type * [n_matrix.out_indexes.size()];
	for (int i = 0; i < new_out_indexes.size(); i++) {
		new_matrix[i] = new Type[n_matrix.in_indexes.size()];
		for (int j = 0; j < new_in_indexes.size(); j++) {
			new_matrix[i][j] = n_matrix.matrix[transformIndex(new_out_indexes[i], n_matrix.out_indexes)][transformIndex(new_in_indexes[j], n_matrix.in_indexes)];
		}
	}
	n_matrix.out_indexes = new_out_indexes;
	n_matrix.in_indexes = new_in_indexes;
	//printf("%i %i", transformIndex(edge[1], n_matrix.out_indexes), transformIndex(edge[1], n_matrix.in_indexes));
	new_matrix[transformIndex(edge[1], n_matrix.out_indexes)][transformIndex(edge[0], n_matrix.in_indexes)] = numeric_limits<Type>::max();
	n_matrix.matrix = new_matrix;
	return n_matrix;
}

template <class Type>
numbered_matrix<Type> clone_numbered_matrix(numbered_matrix<Type> n_matrix) {
	numbered_matrix<Type> new_n_matrix;
	new_n_matrix.matrix = clone_matrix(n_matrix.matrix, n_matrix.out_indexes.size(), n_matrix.in_indexes.size());
	new_n_matrix.out_indexes = vector<Type>(n_matrix.out_indexes);
	new_n_matrix.in_indexes = vector<Type>(n_matrix.in_indexes);
	return new_n_matrix;
}

template<class Type>
void print_numbered_matrix(numbered_matrix<Type> n_matrix) {
	cout << "\t\t";
	for (int i = 0; i < n_matrix.in_indexes.size(); i++) {
		cout << n_matrix.in_indexes.at(i);
		cout << "\t\t";
	}
	cout << "\n";
	for (int i = 0; i < n_matrix.out_indexes.size(); i++) {
		cout << n_matrix.out_indexes.at(i);
		cout << "\t\t";
		for (int j = 0; j < n_matrix.in_indexes.size(); j++) {
			cout << n_matrix.matrix[i][j];
			if (j != n_matrix.in_indexes.size() - 1) {
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
	const int divider = (int)1 / accuracy;

	//генерация координат точек
	float** point_coord = new float* [n];
	for (int i = 0; i < n; i++) {
		point_coord[i] = new float[2];
		point_coord[i][0] = (float)(rand() % (divider + 1)) / divider;
		point_coord[i][1] = (float)(rand() % (divider + 1)) / divider;
	}

	//вычисление матрицы смежности графа на основе расстояний между точками
	float** adj_matrix = new float* [n];
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
				int* pos = new int[2]{ i,j };
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
Type max_zero_degree_null(numbered_matrix<Type> n_matrix, int* pos) {
	Type max = numeric_limits<Type>::min();
	vector<int*> null_positions = matrix_null_element_positions(n_matrix.matrix, n_matrix.out_indexes.size());
	while (!null_positions.empty()) {
		Type val = zero_degree(n_matrix.matrix, null_positions.back(), n_matrix.out_indexes.size());
		if (val > max) {
			max = val;
			pos[0] = n_matrix.out_indexes[null_positions.back()[0]];
			pos[1] = n_matrix.in_indexes[null_positions.back()[1]];
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
		Type** matrix = new Type * [matrix_lines.size()];
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

vector<int> generateIndexList(int n) {
	vector<int> indexes;
	for (int i = 0; i < n; i++) {
		indexes.push_back(i);
	}
	return indexes;
}

template <class Type>
void print_vector(vector<Type> vector) {
	for (int i = 0; i < vector.size(); i++) {
		printf("%i ",vector.at(i));
	}
	printf("\n");
}

int transformIndex(int source_index, vector<int> index_list) { //приведение индексов после удаления рёбер из матрицы смежности
	int new_index;
	new_index = -1;
	for (int i = 0; i < index_list.size(); i++) {
		if (source_index == index_list.at(i)) {
			new_index = i;
		}
	}
	return new_index;
}

template <class T>
struct node_alg_little {
	numbered_matrix<T> adj_matrix;
	int* edge;
	T val;
};

template <class Type>
struct binary_tree_node {
	Type content;
	binary_tree_node* left;
	binary_tree_node* right;
	binary_tree_node* parent;
};

template <class Type>
binary_tree_node<node_alg_little<Type>>* alg_little_step(binary_tree_node<node_alg_little<Type>>* node) {
	int* edge = new int[2];
	max_zero_degree_null(node->content.adj_matrix, edge);

	binary_tree_node<node_alg_little<int>> left_node;
	left_node.content.adj_matrix = remove_edge(clone_numbered_matrix(node->content.adj_matrix), edge);
	left_node.content.edge = edge;
	print_numbered_matrix(left_node.content.adj_matrix);
	cout << endl << endl;
	left_node.content.val = node->content.val + reduce_matrix(left_node.content.adj_matrix.matrix,
		left_node.content.adj_matrix.out_indexes.size(), left_node.content.adj_matrix.in_indexes.size());
	left_node.parent = node;
	node->left = &left_node;

	cout << left_node.content.edge[0] << ' ' << left_node.content.edge[1] << endl;
	cout << left_node.content.val << endl;
	print_numbered_matrix(left_node.content.adj_matrix);
	cout << endl << endl;
	
	binary_tree_node<node_alg_little<int>> right_node;
	right_node.content.adj_matrix = clone_numbered_matrix(node->content.adj_matrix);
	right_node.content.adj_matrix.matrix[transformIndex(edge[0], right_node.content.adj_matrix.out_indexes)]
		[transformIndex(edge[1], right_node.content.adj_matrix.in_indexes)] = numeric_limits<Type>::max();
	right_node.content.edge = edge;
	print_numbered_matrix(left_node.content.adj_matrix);
	cout << endl << endl;
	right_node.content.val = node->content.val + reduce_matrix(right_node.content.adj_matrix.matrix,
		right_node.content.adj_matrix.out_indexes.size(), right_node.content.adj_matrix.in_indexes.size());
	right_node.parent = node;
	node->right = &right_node;

	cout << right_node.content.edge[0] << ' ' << right_node.content.edge[1] << endl;
	cout << right_node.content.val << endl;
	print_numbered_matrix(right_node.content.adj_matrix);
	cout << endl << endl;

	binary_tree_node<node_alg_little<Type>>* nodes = new binary_tree_node<node_alg_little<Type>>[2];
	nodes[0] = left_node;
	nodes[1] = right_node;
	return nodes;
}

template <class Type>
vector<int> traveling_saleman_problem_solution_alg_little(numbered_matrix<Type> adj_matrix) {
	vector<binary_tree_node<node_alg_little<Type>>*> leaf_list;
	binary_tree_node<node_alg_little<Type>> root;
	root.content.adj_matrix = clone_numbered_matrix(adj_matrix);
	root.content.val = reduce_matrix(root.content.adj_matrix.matrix, root.content.adj_matrix.out_indexes.size(),
		root.content.adj_matrix.in_indexes.size());
	root.parent = NULL;

	Type min_val = numeric_limits<Type>::max();
	binary_tree_node<node_alg_little<Type>>* next_node = &root;
	while (next_node->content.adj_matrix.out_indexes.size() + next_node->content.adj_matrix.in_indexes.size() > 2) {
		binary_tree_node<node_alg_little<Type>>* nodes = alg_little_step(next_node);
		next_node->left = &nodes[0];
		next_node->right = &nodes[1];
		leaf_list.push_back(&nodes[0]);
		leaf_list.push_back(&nodes[1]);
		int min_i = 0;
		for (int i = 0; i < leaf_list.size(); i++) {
			if (leaf_list[i]->content.val < min_val) {
				min_val = leaf_list[i]->content.val;
				next_node = leaf_list[i];
				min_i = i;
			}
		}
		leaf_list.erase(leaf_list.begin() + min_i);
		min_val = numeric_limits<Type>::max();
	}
	vector<int*> edges;
	int* edge = new int[2];
	edge[0] = next_node->content.adj_matrix.out_indexes[0];
	edge[1] = next_node->content.adj_matrix.in_indexes[0];
	edges.push_back(edge);
	while (next_node->parent != NULL) {
		if (next_node->parent->right != next_node) {
			edges.push_back(next_node->content.edge);
		}
		next_node = next_node->parent;
	}
	for (int i = 0; i < edges.size(); i++) {
		cout << edges[i][0] << edges[i][1] << endl;
	}
	vector<int> path;
	int next_vertex = 0;
	path.push_back(next_vertex);
	while (edges.empty() != true) {
		for (int i = 0; i < edges.size(); i++) {
			if (edges[i][0] == next_vertex) {
				next_vertex = edges[i][1];
				edges.erase(edges.begin() + i);
				break;
			};
		}
		path.push_back(next_vertex);
	}
	
	return path;
}

int main()
{
	const int n = 5;
	//float** matrix = generateFlatGraph(n);

	//int** matrix = matrix_from_file<int>("input.txt");
	numbered_matrix<int> n_matrix;
	n_matrix.matrix = matrix_from_file<int>("input.txt");
	n_matrix.out_indexes = generateIndexList(n);
	n_matrix.in_indexes = generateIndexList(n);

	vector<int> path = traveling_saleman_problem_solution_alg_little(n_matrix);
	for (int i = 0; i < path.size(); i++) {
		cout << path[i] << ' ';
	}

	/*binary_tree_node<int> root;
	root.val = h;

	n_matrix = remove_edge<int>(n_matrix, edge);
	int h2 = reduce_matrix(n_matrix.matrix, n_matrix.out_indexes.size(), n_matrix.in_indexes.size());
	cout << h2;*/
}

