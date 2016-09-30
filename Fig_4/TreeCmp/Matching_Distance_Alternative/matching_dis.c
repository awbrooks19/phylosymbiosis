/*Author: Vaibhav Rajan, Yu Lin, LCBB
Compute the matching distance for two trees
Input: two tree files t1 and t2, number of leaves in the tree
Assumptions:
- input tree file contains only 1 line which is the tree string, there is no whitespace before or after the string
- input tree is unrooted, binary
- leaf numbering in input file starts from 1
Notes:
- the output file has the name t1.t2.md
- sample intput tree: (((7:0.08,5:0.07):0.16,2:0.27):0.18,(6:0.04,((4:0.00,3:0.00):0.04,8:0.05):0.00):0.01,1:0.08); number of leaves in the tree = 8;
*/


#include <stdio.h>
#include<stdlib.h>
#include <string.h>
#include<limits.h>
#include<math.h>
#include <assert.h>
#include <list>

#include "hungarian.h"

using namespace std;

#define MAXTREESTRLEN 100000
#define MAXLEN 100
#define MAXLEAF 10000

struct Ptree
{
	int leaf_number;
	int parent[MAXLEAF];
	int lchild[MAXLEAF];
	int rchild[MAXLEAF];
	int edge[MAXLEAF][MAXLEAF]; //vector representation for bipatitions
};



struct Ptree tree1, tree2;

void pttree(struct Ptree *treeA, int node){
	int i;
 	//getchar();
	//int total = treeA.leaf_number*2-1;
//printf("Node:%d, leaveNUM: %d \n", node, (*treeA).leaf_number);getchar();
	if (((*treeA).lchild[node] == -1) &&  ((*treeA).rchild[node] == -1)){
		if (node != 0){
			for(i=0; i<(*treeA).leaf_number; i++){	
				(*treeA).edge[node][i] = 0;
			}
			(*treeA).edge[node][node] = 1;
		}else{
			for(i=0; i<(*treeA).leaf_number; i++){	
				(*treeA).edge[node][i] = 1;
			}
			(*treeA).edge[node][node] = 0;
		}
		//printf("LEAVE Node:%d\n", node);getchar();
	}else{
		//printf("Node:%d -> (%d,%d)\n", node,(*treeA).lchild[node],(*treeA).rchild[node]);getchar();
		if ((*treeA).lchild[node] != -1)
			pttree(treeA, (*treeA).lchild[node]);
		if ((*treeA).rchild[node] != -1)
			pttree(treeA, (*treeA).rchild[node]);
		for(i=0; i<(*treeA).leaf_number; i++){	
			(*treeA).edge[node][i] = ( (*treeA).edge[(*treeA).lchild[node]][i] + (*treeA).edge[(*treeA).rchild[node]][i] ) % 2;
		}
		//printf("INTER Node:%d\n", node);getchar();
	}
}


int vecdistance(struct Ptree *treeA, int a, struct Ptree *treeB, int b){
	int i;	
	for(i=0; i<(*treeA).leaf_number; i++){	
		if ( (*treeA).edge[a][i] != (*treeB).edge[b][i] )
			return(0);
	}
	return(1);
}

int rfdistance(struct Ptree *treeA, struct Ptree *treeB){
	int i,j;
	int rf = 0;
	int tag = 0;
	int flag[12000];
	int ii;
	double sum;
	sum = 0;
	for(j=(*treeB).leaf_number+1; j<(*treeB).leaf_number * 2-1; j++){
		flag[j] = 0;
	}

	for(i=(*treeA).leaf_number+2; i<(*treeA).leaf_number * 2-1; i++){	
		tag = 0;
		//printf("%d\n",i);
		for(j=(*treeB).leaf_number+1; ((j<(*treeB).leaf_number * 2-1)&&(!tag)); j++){
			if (vecdistance(treeA, i, treeB, j) > 0){
				tag++;
				rf++;
				flag[j]++;
			}
		}
	}
	rf = (*treeA).leaf_number - 3 - rf;
	return(rf);
}


void compute_matrix(int *r, int range){
       int e1,e2;
       int i;
       int temp1;
       int temp2;
       int row_r, col_r;
       row_r = 0;
       col_r = 0;
       int count = 0;
       int spe1;
       int spe2;
       if (tree1.lchild[tree1.leaf_number] > tree1.leaf_number)
               spe1 = tree1.lchild[tree1.leaf_number];
       else
               spe1 = tree1.rchild[tree1.leaf_number];
       if (tree2.lchild[tree2.leaf_number] > tree2.leaf_number)
               spe2 = tree2.lchild[tree2.leaf_number];
       else
               spe2 = tree2.rchild[tree2.leaf_number];
       for(e1 = tree1.leaf_number + 1; e1 < 2*tree1.leaf_number - 1; e1++){
               if (e1 != spe1){
                       col_r = 0;
                       for(e2 = tree2.leaf_number + 1; e2 < 2*tree2.leaf_number - 1; e2++){
                               if (e2 != spe2){
                                       temp1 = 0;
                                       temp2 = 0;
                                       for(i=0; i<tree2.leaf_number; i++){
                                               if ( tree1.edge[e1][i] != tree2.edge[e2][i] )
                                                       temp1++;
                                       }
                                       for(i=0; i<tree2.leaf_number; i++){
                                               if ( tree1.edge[e1][i] == tree2.edge[e2][i] )
                                                       temp2++;
                                       }
                                       if (temp1 < temp2){
                                                       r[row_r*range+col_r] = temp1;
                                       }else{
                                                       r[row_r*range+col_r] = temp2;
                                       }
                                       col_r++;
                               }
                       }
			if ( col_r !=  tree1.leaf_number - 3 ){
				printf("overflow of col_r! %d\n", col_r); getchar();
			}
                       row_r++;
               }
       }
	if (row_r !=  tree1.leaf_number - 3) {
		printf("overflow of row_r! %d\n", row_r); getchar();
	}
}

int** array_to_matrix(int* m, int rows, int cols){
  int i,j;
  int** r;
  r = (int**)calloc(rows,sizeof(int*));
  if (r == NULL){
	printf("alloc r failed!\n"); getchar();
  }
  for(i=0;i<rows;i++)
  {
    r[i] = (int*)calloc(cols,sizeof(int));
	  if (r[i] == NULL){
		printf("alloc ri failed!\n"); getchar();
	  }
    for(j=0;j<cols;j++)
      r[i][j] = m[i*cols+j];
  }
  return r;
}


void newick2lcbb(char *input_filename, int num_leaves, struct Ptree *tree)
{
	FILE *fp, *fq;
	char tree_str[MAXTREESTRLEN];
	int root;

	fp = fopen(input_filename, "r");
	fgets(tree_str, MAXTREESTRLEN, fp);	
	fclose(fp);

	//printf("%d %s\n", strlen(tree_str), tree_str);

	list<int> l;
	list<double> ll;
	char temp[MAXLEN];
	int state = 0;
	int index;
	int int_node = num_leaves + 2;
	int node1, node2;
	double edgelen;
	int i,j;
	
	root = num_leaves+1;
	(*tree).leaf_number = num_leaves;
	for(i=0; i<num_leaves*2; i++){	
		(*tree).parent[i] = -1;
		(*tree).lchild[i] = -1;
		(*tree).rchild[i] = -1;
		for(j=0; j<num_leaves*2; j++){
			(*tree).edge[i][j] = 0;
		}
	}

	for(i = 0; i < strlen(tree_str); i++){
		if (tree_str[i] == ')'){
			if (state == 3){
				int k = 0;
				for(int j = index; j < i; j++, k++)	
					temp[k] = tree_str[j];
				temp[k] = 0;
				ll.push_front(atof(temp));
				//printf(") pushing %lf\n", atof(temp));
				//getchar();
				state = 3;
			}
			node1 = l.front();
			l.pop_front();
			edgelen = ll.front(); 
			ll.pop_front();
			//edge
			//printf("%d %d %lf\n", int_node, node1, edgelen);
			(*tree).parent[node1-1] = int_node-1;
			if ((*tree).lchild[int_node-1] == -1){
				//tree2.lchild[parent-1] = child-1;
				(*tree).lchild[int_node-1] = node1-1;
			}else{ 
				//tree2.rchild[parent-1] = child-1;
				(*tree).rchild[int_node-1] = node1-1;
			}

			node1 = l.front();
			l.pop_front();
			edgelen = ll.front();
			ll.pop_front();
			//edge
			//printf("%d %d %lf\n", int_node, node1, edgelen);
			(*tree).parent[node1-1] = int_node-1;
			if ((*tree).lchild[int_node-1] == -1){
				//tree2.lchild[parent-1] = child-1;
				(*tree).lchild[int_node-1] = node1-1;
			}else{ 
				//tree2.rchild[parent-1] = child-1;
				(*tree).rchild[int_node-1] = node1-1;
			}
			l.push_front(int_node);
			//printf(") pushing %d\n", int_node);
			//getchar();
			int_node += 1;	
		}
		else if (tree_str[i] == ':'){
			if (state == 1){	
				int k = 0;	
				for(int j = index; j < i; j++, k++)
                                        temp[k] = tree_str[j];
                                temp[k] = 0;
				l.push_front(atoi(temp));
				//printf(": pushing %d\n", atoi(temp)+1);
				//getchar();
			}
			state = 2;
		}
		else if (tree_str[i] == ','){
			if (state == 3){
				int k = 0;
				for(int j = index; j < i; j++, k++)
                                        temp[k] = tree_str[j];
                                temp[k] = 0;
				ll.push_front(atof(temp));
				//printf(", pushing %lf\n", atof(temp));
				//getchar();
				state = 0;
			}	
		}
		else if (tree_str[i] != '('){
			if (state == 0){
				index = i;
				state = 1;
			}
			else if (state == 2){
				index = i;
				state = 3;
			}
		}
	}
	if (l.size() == 2){
		node1 = l.front();
		l.pop_front();
		node2 = l.front();
		l.pop_front();
		//edgelen = ll.front();
		//ll.pop_front();
		//fake root edges
		//printf("%d %d\n", num_leaves+1, node1);
		//printf("%d %d\n", num_leaves+1, node2);
		(*tree).parent[node1-1] = num_leaves;
		if ((*tree).lchild[num_leaves] == -1){
			(*tree).lchild[num_leaves] = node1-1;
		}else{ 
			(*tree).rchild[int_node-1] = node1-1;
		}
		(*tree).parent[node2-1] = num_leaves;
		if ((*tree).lchild[num_leaves] == -1){
			(*tree).lchild[num_leaves] = node2-1;
		}else{ 
			(*tree).rchild[num_leaves] = node2-1;
		}
	}
	else{
		printf("input format incorrect! compare with sample tree of 8 leaves: (((6:0.08,4:0.07):0.16,1:0.27):0.18,(5:0.04,((3:0.00,2:0.00):0.04,7:0.05):0.00):0.01,0:0.08);\n");getchar();
	}
	//printf("input correct!\n"); getchar();
	pttree(tree,num_leaves);

}


int trees_mmdis(int num_leaf){
	int root, root2;
	int i,j;
	int child, parent;
	int* r;
	int** m;
	int mmdis;
	int step;
	int tempstep;
	double dtemp;
	hungarian_problem_t p;

	root = num_leaf;
	pttree(&tree1,root);
	pttree(&tree2,root);	
	r = (int *)malloc(sizeof(int)*((root - 3)*(root - 3)));
	if (r == NULL){
		printf("malloc r failed!\n"); getchar();
	}
	compute_matrix(r,root-3);
	m = array_to_matrix(r,root-3,root-3);
	  /* initialize the gungarian_problem using the cost matrix*/
	hungarian_init(&p, m , root-3, root-3, HUNGARIAN_MODE_MINIMIZE_COST) ;
  	/* solve the assignement problem */
	mmdis = hungarian_solve(&p);
	  /* free used memory */
	hungarian_free(&p);
	for(i=0;i<root-3;i++){
		free(m[i]);
	}
	free(m);
	free(r);
	return(mmdis);
}



int main(int argc, char *argv[]){
        FILE *f1, *f2,*f_rf,*f_mm;
	char file_rf[100], file_mm[100];
	int r1,r2;
	int i;
	int nol;
	if(argc != 4){
		printf("Usage: matching_dis treefile1 treefile2 num_of_leaves\n");
		exit(1);
	}
	nol = atoi(argv[3]);	
	sprintf(file_mm,  "%s.%s_md", argv[1], argv[2]);
	f_mm = fopen(file_mm, "w");
	newick2lcbb(argv[1], nol, &tree1);
	//printf("tree1 complete\n"); getchar();
	newick2lcbb(argv[2], nol, &tree2);
	//printf("tree2 complete!\n"); getchar();
	fprintf(f_mm, "%d\n", trees_mmdis(nol)); 
	//printf("%d\n", trees_mmdis(nol)); 
	fclose(f_mm);

}
