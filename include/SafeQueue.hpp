// https://wangpengcheng.github.io/2019/05/17/cplusplus_theadpool/
#pragma once 
#include <mutex> 
#include <queue> 
// Thread safe implementation of a Queue using a std::queue

template <typename T>
class SafeQueue {
private:
  std::queue<T> m_queue; // Constructing a queue using template functions

  std::mutex m_mutex;  // Accessing a mutex semaphore

public:
  SafeQueue() { // Empty Constructor


  }

  SafeQueue(SafeQueue& other) {// Copy Constructor

    //TODO:
  }

  ~SafeQueue() { // Destructor

  }


  bool empty() {  // Is the queue empty?

    std::unique_lock<std::mutex> lock(m_mutex); // The mutex signal variable is locked to prevent m_queue from being changed

    return m_queue.empty();
  }
  
  int size() {
    std::unique_lock<std::mutex> lock(m_mutex); // The mutex signal variable is locked to prevent m_queue from being changed

    return m_queue.size();
  }

  // Adding elements to a queue
  void enqueue(T& t) {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_queue.push(t);
  }

  // Remove elements from the queue
  bool dequeue(T& t) {
    std::unique_lock<std::mutex> lock(m_mutex); // Queue Lock

    if (m_queue.empty()) {
      return false;
    }
    t = std::move(m_queue.front()); // Take out the first element of the team, return the value of the first element of the team, and make a right value reference
    
    m_queue.pop(); // Pop the first element into the queue

    return true;
  }
};