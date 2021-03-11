#include <iostream>
#include <math.h>
#include <limits>
#include <vector>
#include <type_traits>
#include <fstream>
#include <string>
using namespace std;

const float accuracy = 0.0001;

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
	if (transformIndex(edge[1], n_matrix.out_indexes) != -1 && transformIndex(edge[0], n_matrix.in_indexes) != -1) {
		new_matrix[transformIndex(edge[1], n_matrix.out_indexes)][transformIndex(edge[0], n_matrix.in_indexes)] = numeric_limits<Type>::max();
	}
	n_matrix.matrix = new_matrix;
	return n_matrix;
}

template <class Type>
numbered_matrix<Type> clone_numbered_matrix(numbered_matrix<Type> n_matrix) {
	numbered_matrix<Type> new_n_matrix;
	new_n_matrix.matrix = clone_matrix(n_matrix.matrix, n_matrix.out_indexes.size(), n_matrix.in_indexes.size());
	new_n_matrix.out_indexes = vector<int>(n_matrix.out_indexes);
	new_n_matrix.in_indexes = vector<int>(n_matrix.in_indexes);
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
	const int divider = (int)1 / sqrt(accuracy);

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
	for (int i = 0; i < null_positions.size(); i++){
		Type val = zero_degree(n_matrix.matrix, null_positions[i], n_matrix.out_indexes.size());
		if (val > max) {
			max = val;
			pos[0] = n_matrix.out_indexes[null_positions[i][0]];
			pos[1] = n_matrix.in_indexes[null_positions[i][1]];
		}
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
string string_node_state(binary_tree_node<node_alg_little<Type>>* node) {
	string str = "";
	if (node->left != NULL) {
		str += "[ " + string_node_state(node->left) + " ] <";
	}
	str += " ";
	if (node->parent == NULL) {
		str += "* ";
	}
	if (node->content.edge != NULL) {
		str += to_string(node->content.edge[0]) + "-" + to_string(node->content.edge[1]) + " ";
	}
	str += "(" + to_string(node->content.val) + ")";
	if (node->parent == NULL) {
		str += " *";
	}
	if (node->right != NULL) {
		str += " > [ " + string_node_state(node->right) + " ]";
	}
	return str;
}

template <class Type>
vector<int*> get_edges(binary_tree_node<node_alg_little<Type>>* node) {
	vector<int*> edges;
	binary_tree_node<node_alg_little<Type>>* curr_node = node;
	while (curr_node->parent != NULL) {
		if (curr_node->parent->right != curr_node) {
			edges.push_back(curr_node->content.edge);
		}
		curr_node = curr_node->parent;
	}
	return edges;
}

vector<int> get_cycle(vector<int*> edges) {
	vector<int> path;
	int start_vertex = edges[0][0];
	int next_vertex = edges[0][1];
	path.push_back(start_vertex);
	while (next_vertex != start_vertex && !(path.size() > edges.size() + 2)) {
		for (int i = 0; i < edges.size(); i++) {
			if (edges[i][0] == next_vertex) {
				path.push_back(next_vertex);
				next_vertex = edges[i][1];
				break;
			};
			if (i == edges.size() - 1) {
				vector<int> no_cycle;
				no_cycle.push_back(-1);
				return no_cycle;
			}
		}
	}
	path.push_back(next_vertex);

	/*cout << "found cycle: ";
	for (int i = 0; i < path.size(); i++) {
		cout << path[i] << " ";
	}
	cout << endl;*/

	return path;
}

template <class Type>
binary_tree_node<node_alg_little<Type>>* alg_little_step(binary_tree_node<node_alg_little<Type>>* node) {
	int* edge = new int[2];
	max_zero_degree_null(node->content.adj_matrix, edge);
	//cout << edge[0] << " " << edge[1] << " ";

	binary_tree_node<node_alg_little<Type>> left_node;
	left_node.content.adj_matrix = remove_edge(clone_numbered_matrix(node->content.adj_matrix), edge);
	left_node.content.edge = edge;
	left_node.content.val = node->content.val + reduce_matrix(left_node.content.adj_matrix.matrix,
		left_node.content.adj_matrix.out_indexes.size(), left_node.content.adj_matrix.in_indexes.size());
	//cout << "in: " << left_node.content.val << " ";
	left_node.left = NULL;
	left_node.right = NULL;
	left_node.parent = node;
	node->left = &left_node;

	/*cout << "+ " << left_node.content.edge[0] << ' ' << left_node.content.edge[1] << endl;
	cout << left_node.content.val << endl;
	print_numbered_matrix(left_node.content.adj_matrix);
	cout << endl << endl;*/
	
	binary_tree_node<node_alg_little<Type>> right_node;
	right_node.content.adj_matrix = clone_numbered_matrix(node->content.adj_matrix);
	right_node.content.adj_matrix.matrix[transformIndex(edge[0], right_node.content.adj_matrix.out_indexes)]
		[transformIndex(edge[1], right_node.content.adj_matrix.in_indexes)] = numeric_limits<Type>::max();
	right_node.content.edge = edge;
	right_node.content.val = node->content.val + reduce_matrix(right_node.content.adj_matrix.matrix,
		right_node.content.adj_matrix.out_indexes.size(), right_node.content.adj_matrix.in_indexes.size());
	//cout << "ex: " << right_node.content.val << endl;
	right_node.left = NULL;
	right_node.right = NULL;
	right_node.parent = node;
	node->right = &right_node;

	/*cout << "- " << right_node.content.edge[0] << ' ' << right_node.content.edge[1] << endl;
	cout << right_node.content.val << endl;
	print_numbered_matrix(right_node.content.adj_matrix);
	cout << endl << endl;*/

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
	root.content.edge = NULL;
	root.parent = NULL;

	Type min_val = numeric_limits<Type>::max();
	binary_tree_node<node_alg_little<Type>>* next_node = &root;
	while (next_node->content.adj_matrix.out_indexes.size() + next_node->content.adj_matrix.in_indexes.size() > 2) {
		binary_tree_node<node_alg_little<Type>>* nodes = alg_little_step(next_node);
		next_node->left = &nodes[0];
		next_node->right = &nodes[1];
		//cout << string_node_state(&root) << endl;
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
		/*if (next_node->parent->left == next_node) {
			cout << "selected in: ";
		}
		else {
			cout << "selected ex: ";
		}
		cout << next_node->content.edge[0] << " " << next_node->content.edge[1] << " val: " << next_node->content.val << endl;*/
		//исключение рёбер которые могут состовлять подциклы с текщими рёбрами
		if (next_node->parent->left == next_node) {
			vector<int*> edges = get_edges(next_node);
			for (int i = 0; i < edges.size(); i++) {
				for (int j = 0; j < edges.size(); j++) {
					if (i != j) {
						vector<int*> edges_cycle(edges);
						int* test_edge = new int[2];
						test_edge[0] = edges[i][1];
						test_edge[1] = edges[j][0];
						edges_cycle.push_back(test_edge);
						vector<int> test_result = get_cycle(edges_cycle);
						if (test_result[0] != -1 && test_result.size() != root.content.adj_matrix.out_indexes.size() + 1) {
							//cout << "blocked:" << test_edge[0] << " " << test_edge[1] << endl;
							if (transformIndex(edges[i][1], next_node->content.adj_matrix.out_indexes) != -1 &&
								transformIndex(edges[j][0], next_node->content.adj_matrix.in_indexes) != -1) {
								next_node->content.adj_matrix.matrix[transformIndex(edges[i][1], next_node->content.adj_matrix.out_indexes)]
									[transformIndex(edges[j][0], next_node->content.adj_matrix.in_indexes)] = numeric_limits<Type>::max();
							}
						}
					}
				}
			}
		}
		leaf_list.erase(leaf_list.begin() + min_i);
		min_val = numeric_limits<Type>::max();
	}

	vector<int*> edges = get_edges(next_node);
	int* edge = new int[2];
	edge[0] = next_node->content.adj_matrix.out_indexes[0];
	edge[1] = next_node->content.adj_matrix.in_indexes[0];
	edges.push_back(edge);
	
	/*for (int i = 0; i < edges.size(); i++) {
		cout << edges[i][0] << edges[i][1] << endl;
	}*/

	return get_cycle(edges);
}

template <class Type>
Type path_len(vector<int> path, Type** matrix) {
	Type sum = 0;
	for (int i = 0; i < path.size() - 1; i++) {
		sum += matrix[path[i]][path[i + 1]];
	}
	return sum;
}

int main()
{

	//матрица из файла input.txt
	const int n = 5;
	numbered_matrix<int> n_matrix;
	n_matrix.matrix = matrix_from_file<int>("input2.txt");

	//случайно сгенерированная матрица
	/*const int n = 20;
	numbered_matrix<float> n_matrix;
	srand(0);
	n_matrix.matrix = generateFlatGraph(n);*/

	n_matrix.out_indexes = generateIndexList(n);
	n_matrix.in_indexes = generateIndexList(n);

	print_numbered_matrix(n_matrix);

	vector<int> path = traveling_saleman_problem_solution_alg_little(n_matrix);

	float val = path_len(path, n_matrix.matrix);

	for (int i = 0; i < path.size(); i++) {
		cout << path[i] << ' ';
	}
	cout << endl;
	cout << "val: " << val;
}

