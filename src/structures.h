#ifndef AFFINE_STRUCTURES_H
#define AFFINE_STRUCTURES_H

#include <stdbool.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * In structures, you will find all that is needed/used for the different structures.
 */

/**
 * A truth table which holds the information of the n of the truth table and all its elements.
 */
typedef struct TruthTable {
    size_t n; // Dimension of the function
    size_t *elements; // All the (2^n) elements of the function
} TruthTable;

/**
 * Initialize a new TruthTable. This method takes in the dimension (n) of the function and creates a new TruthTable where
 * the list for the elements is malloced with 2^n, where n is the n of the function.
 * The dimension and a empty list for the elements is declared.
 * @param n The n of the function
 * @return The pointer to the new TruthTable
 */
TruthTable *initTruthTable(size_t n);

/**
 * Add the elements from one function into another.
 * @param dest The truth table to add to
 * @param src The truth table with the elements to add to the other.
 */
void add(TruthTable *dest, TruthTable *src);

/**
 * Compose two functions together and return a new function containing the result
 * @param f The function F
 * @param g The function G
 * @return H = F * G
 */
TruthTable *compose(TruthTable *f, TruthTable *g);

/**
 * Find the inverse of F, F^{-1}
 * @param f The function F
 * @return F inverse
 */
TruthTable * inverse(TruthTable *f);

/**
 * Create a random affine function with 2^n elements
 * @param n The dimension
 * @return A new random Affine function
 */
TruthTable *randomLinearFunction(size_t n);

/**
 * Create a random affine permutation with 2^n elements
 * @param n The dimension
 * @return A new random Affine permutation
 */
TruthTable *randomLinearPermutation(size_t n);

/**
 * Create a new affine truth table with respect to a function F
 * @param F The function F
 * @return A new affine function G
 */
TruthTable *createAffineTruthTable(TruthTable *F);

/**
 * Create new linear truth table with respect to a function F where F and G are EA-equivalent
 * @param F The function F
 * @return A new linear EA function G
 */
TruthTable *createLinearFunction(TruthTable *F);

/**
 * Print all the elements of the TruthTable to the console
 * @param tt The pointer to the truth table to print out
 */
void printTruthTable(TruthTable *tt);

/**
 * Free the memory that is allocated for the struct TruthTable
 * @param tt The pointer for the TruthTable that is to be destroyed
 */
void destroyTruthTable(TruthTable *tt);

/**
 * Holds all the information for the Partition:
 * Number of buckets, all the multiplicities, the sizes of the buckets, and all the elements in all the buckets.
 *
 */
typedef struct Partition {
    size_t numBuckets; // Number of buckets
    size_t *multiplicities; // All the multiplicities that matches the buckets, this is a list of lists
    size_t *bucketSizes; // A list that holds the information of the sizes of all the buckets
    size_t **buckets; // A list of lists, holds the elements in each bucket.
} Partition;

/**
 * Initialize a new Partition. This method takes in the dimension (n) of the function and created a new Partition where
 * it allocates all the memory needed to create all the lists.
 * @param n // The dimension of the field we are working in
 * @return A pointer to a new Partition
 */
Partition *initPartition(size_t n);

/**
 * Print out the elements to each bucket to the console. Each line represents a bucket.
 * @param partition The partition to print out
 */
void printPartitionBuckets(Partition *partition);

/**
 * Print all the pre images of the partitions
 * @param F The partition to print the images of
 */
void printPartition(Partition *F);

/**
 * Perform the partitioning of a function. This function will find out which bucket partition each element from the
 * TruthTable belongs to.
 * @param tt The truth table to partition
 * @return A new Partition of the truth table.
 */
Partition *partitionTt(TruthTable *tt);

/**
 * Free the memory allocated for the Partition
 * @param partition The Partition to destroy
 */
void destroyPartition(Partition *partition);

/**
 * A linked list for holding numbers as the data
 */
typedef struct Node {
    size_t data; // The data for this node
    struct Node *next; // Pointer to the next node in the linked list
} Node;

/**
 * Initial a new empty Linked list for numbers.
 * Note: The HEAD of the list is a empty node, and will always be empty.
 * The new nodes is inserted in the head.next, resulting in that the list is in reverse.
 * Ex: Suppose you insert 3 new nodes (1,2,3) to a empty list, the list will look like: [head -> 3 -> 2 -> 1]
 * @return A empty head node for a linked list
 */
Node *initNode();

/**
 * Add a new node to a linked list
 * Note: The HEAD of the list is a empty node, and will always be empty.
 * The new nodes is inserted in the head.next, resulting in that the list is in reverse.
 * Ex: Suppose you insert 3 new nodes (1,2,3) to a empty list, the list will look like: [head -> 3 -> 2 -> 1]
 * @param head
 * @param data
 */
void addNode(Node *head, size_t data);

/**
 * Get the data from the node requested.
 * The linked list is in reverse and the head node is a empty node, s.t. the list is on the form [head -> 3 -> 2 -> 1]
 * where the size of the linked list is 3
 * @param head The pointer to the head of the list
 * @param index The element to get
 * @return The number for the element in position = index
 */
size_t getNode(Node *head, size_t index);

/**
 * Get the size of the linked list
 * @param head The pointer to the head of the list
 * @return The size of the list
 */
size_t countNodes(Node *head);

/**
 * Print out all the nodes in a linked list
 * @param head The pointer to the head of the linked list
 */
void printNodes(Node *head);

/**
 * Free all allocated memory for a linked list Node
 * @param head The pointer to the head of the linked list
 */
void destroyNodes(Node *head);

/**
 * This is a linked list of TruthTables
 */
typedef struct TtNode {
    TruthTable *data;
    struct TtNode *next; // Pointer to the next node
} TtNode;

/**
 * Initial a new empty Linked list of TruthTables
 * Note: The HEAD of the list is a empty node, and will always be empty.
 * The new nodes is inserted in the head.next, resulting in that the list is in reverse.
 * Ex: Suppose you insert 3 new nodes (1,2,3) to a empty list, the list will look like: [head -> 3 -> 2 -> 1]
 * @return A empty head node for a linked list
 */
TtNode *initTtNode();

/**
 * Add a new TruthTable to the linked list
 * Note: The HEAD of the list is a empty node, and will always be empty.
 * The new nodes is inserted in the head.next, resulting in that the list is in reverse.
 * Ex: Suppose you insert 3 new nodes (1,2,3) to a empty list, the list will look like: [head -> 3 -> 2 -> 1]
 * @param head The head of the linked list
 * @param data The truth table to add
 */
void addTtNode(TtNode *head, TruthTable *data);

/**
 * Get the data from the node requested.
 * The linked list is in reverse, s.t. the list is on the form [1 -> 3 -> 2]
 * where the size of the linked list is 3
 * @param head The head of the linked list
 * @param index The position in the linked list to get the truth table
 * @return The truth table in the position index in the linked list
 */
TruthTable *getTtNode(TtNode *head, size_t index);

/**
 * Get the size of the linked list
 * @param head The pointer to the head of the list
 * @return The size of the list
 */
size_t countTtNodes(TtNode *head);

/**
 * Free all allocated memory for a linked list
 * @param head The pointer to the head of the linked list
 */
void destroyTtNode(TtNode *head);

/**
 * A struct for keeping all the run times for the entire program
 */
typedef struct RunTimes {
    double total;
} RunTimes;

/**
 * Initialize a new empty struct of Run times
 * @return A empty RunTimes
 */
RunTimes *initRunTimes();

/**
 * Stop the time and return the run time
 * @param runTime The time to do the calculations on
 * @param startParsing The start time
 * @return The total run time
 */
double stopTime(double runTime, clock_t startParsing);

/**
 * Print all times to the console
 * @param runTimes the Run times to print
 */
void printTimes(RunTimes *runTimes);

/**
 * Destroy the run times
 * @param runTimes The run times to destroy
 */
void destroyRunTimes(RunTimes *runTimes);

/**
 * Create a standard basis, {b_1, ..., b_n}
 * @param n The dimension
 * @return Standard basis
 */
size_t *createStandardBasis(size_t n);

/**
 * Add a constant c to a function F, s.t. F = F + c
 * @param F A function F to add a constant to
 * @param c The value of the constant
 */
void addConstant(TruthTable *F, size_t c);

#endif //AFFINE_STRUCTURES_H
