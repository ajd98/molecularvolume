#include "queue.h"

// encapsulate grid indices
struct Triple
{
    int x, y, z;
};

// fifo data structure
struct Queue 
{
    int front, rear, size;
    unsigned int capacity;
    Triple* items;
}

// Make a new queue with maximum capacity of ``capacity`` Triples
int newQueue(int capacity)
{
    Queue* queue;
    queue->front = 0;
    queue->size = 0;
    queue->rear = 0;
    queue->capacity = capacity;
    queue->items = (Triple*)malloc(capacity*sizeof(Triple));
}

// Check if the queue is full
int isFull(Queue* queue)
{
    if (queue->capacity == queue->size){
        return 1;
    } else {
        return 0;
    }
}

// Add a Triple to the queue
int enqueue(Triple coords, Queue* queue)
{
    if isFull(queue){
        return 1; // error; queue is full
    }
    // Add the item to the rear of the queue
    // wrapping around to the beginning of the array if need be
    queue->items[queue->rear % queue->capacity]  = coords;
    queue->rear = (queue->rear+1) % queue->capacity;
    queue->size += 1;
    return 0;
}

int isEmpty(Queue* queue)
{
    if (queue->size == 0){
        return 1;
    } else {
        return 0;
    }
    
}

Triple dequeue(Queue* queue)
{
    if isEmpty(queue){
        return NULL; // queue is empty; return a null pointer
    }
    Triple item = queue->items[queue->front];
    queue->front = (queue->front+1)%queue->capacity;
    queue->size -= 1;
    return item;
}

int main()
{
    queue = newQueue(10);
    Triple item1;
    item1.x = 1
    item1.y = 1
    item1.z = 1
    enqueue(item1, queue)
}
