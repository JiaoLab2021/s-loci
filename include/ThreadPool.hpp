// https://wangpengcheng.github.io/2019/05/17/cplusplus_theadpool/
#pragma once 
#include <functional> 
#include <future> 
#include <mutex> 
#include <queue> 
#include <thread> 
#include <utility> 
#include <vector> 
#include "SafeQueue.hpp" 

class ThreadPool {
private:
  class ThreadWorker {  // Built-in thread worker class

  private:
    int m_id;  // Work ID

    ThreadPool * m_pool;  // Thread pool

  public:
    // Constructor

    ThreadWorker(ThreadPool * pool, const int id) 
      : m_pool(pool), m_id(id) {
    }
    // Overriding the `()` operation

    void operator()() {
      std::function<void()> func;  // Define the base function class func
      
      bool dequeued;  // Is the element being taken out of the queue?
      
      // Determine whether the thread pool is closed. If not, loop to extract

      while (!m_pool->m_shutdown) {
        {
          // Lock the thread environment lock to access the sleep and wake-up of the worker thread

          std::unique_lock<std::mutex> lock(m_pool->m_conditional_mutex);
          // If the task queue is empty, block the current thread

          if (m_pool->m_queue.empty()) {
            m_pool->m_conditional_lock.wait(lock); // Wait for condition variable notification and start thread

          }
          // Remove the element from the task queue

          dequeued = m_pool->m_queue.dequeue(func);
        }
        // If the extraction is successful, execute the work function

        if (dequeued) {
          func();
        }
      }
    }
  };

  bool m_shutdown;  // Is the thread pool closed?

  SafeQueue<std::function<void()>> m_queue;  // Execute function security queue, that is, task queue

  std::vector<std::thread> m_threads;  // Worker thread queue

  std::mutex m_conditional_mutex;  // Thread sleep lock mutex variable

  std::condition_variable m_conditional_lock;  // Thread environment lock allows threads to sleep or wake up

public:
    // Thread Pool Constructor

  ThreadPool(const int n_threads)
    : m_threads(std::vector<std::thread>(n_threads)), m_shutdown(false) {
  }

  ThreadPool(const ThreadPool &) = delete;  // Copy the constructor and cancel the default parent class constructor

  ThreadPool(ThreadPool &&) = delete;  // Copy constructor, allowing rvalue references

  ThreadPool & operator=(const ThreadPool &) = delete;  // Assignment Operation

  ThreadPool & operator=(ThreadPool &&) = delete;  // Assignment Operation

  // Inits thread pool

  void init() {
    for (int i = 0; i < m_threads.size(); ++i) {
      m_threads[i] = std::thread(ThreadWorker(this, i));  // Allocating Worker Threads

    }
  }

  // Waits until threads finish their current task and shutdowns the pool

  void shutdown() {
    m_shutdown = true;
    m_conditional_lock.notify_all();  // Notify wake up all worker threads
    
    for (int i = 0; i < m_threads.size(); ++i) {
      if(m_threads[i].joinable()) {  // Determine whether the thread is waiting

        m_threads[i].join();   // Add the thread to the waiting queue

      }
    }
  }

  // Submit a function to be executed asynchronously by the pool
  // The main working function of the thread uses the post-return type to automatically determine the function return value

  template<typename F, typename...Args>
  auto submit(F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
    // Create a function with bounded parameters ready to execute
    // 

    std::function<decltype(f(args...))()> func = std::bind(std::forward<F>(f), std::forward<Args>(args)...);  // Connect function and parameter definitions, special function types, avoid left and right value errors

    // Encapsulate it into a shared ptr in order to be able to copy construct // assign 
    // Encapsulate the task object to facilitate another thread to view the result

    auto task_ptr = std::make_shared<std::packaged_task<decltype(f(args...))()>>(func);

    // Wrap packaged task into void function
    // Using a regular expression, return a function object

    std::function<void()> wrapper_func = [task_ptr]() {
      (*task_ptr)(); 
    };

    // Queue universal security package function and push it into the security queue

    m_queue.enqueue(wrapper_func);

    // Wake up a waiting thread

    m_conditional_lock.notify_one();

    // Returns a previously registered task pointer

    return task_ptr->get_future();
  }

  // Returns the number of queued tasks
  auto get_queue() {
    return m_queue.size();
  }
};