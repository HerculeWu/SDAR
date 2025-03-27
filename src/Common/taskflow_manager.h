#pragma once

#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/for_each.hpp>
#include <thread>
#include <memory>
#include <mutex>
#include <iostream>
namespace TF
{
    class Manager
    {
    public:
        static tf::Taskflow &get_taskflow()
        {
            std::call_once(init_flag_, []()
                           { taskflow_ = std::make_unique<tf::Taskflow>(); });
            return *taskflow_;
        }

        static tf::Executor &get_executor(size_t num_threads = 0)
        {
            std::call_once(exec_init_flag_, [num_threads]()
                           {
            size_t threads_to_use = num_threads;
    
            if (threads_to_use == 0) {
            threads_to_use = std::thread::hardware_concurrency();
            if (threads_to_use == 0) {
                threads_to_use = 4; // fallback 默认值
                std::cerr << "[Taskflow][Warning] Could not detect hardware concurrency. Falling back to 4 threads.\n";
            }
            }
    
            executor_ = std::make_unique<tf::Executor>(256*1024*1024, threads_to_use); });
            return *executor_;
        }

        Manager(const Manager &) = delete;
        Manager &operator=(const Manager &) = delete;

    private:
        Manager() = default;

        static std::once_flag init_flag_;
        static std::once_flag exec_init_flag_;
        static std::unique_ptr<tf::Taskflow> taskflow_;
        static std::unique_ptr<tf::Executor> executor_;
    };

    inline std::once_flag Manager::init_flag_;
    inline std::once_flag Manager::exec_init_flag_;
    inline std::unique_ptr<tf::Taskflow> Manager::taskflow_;
    inline std::unique_ptr<tf::Executor> Manager::executor_;
}