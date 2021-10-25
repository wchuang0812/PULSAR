# Load balance strategy
# Input : Command_runtime and cores_number.
# Command_runtime is a list storing each command's runtime in single thread. 
# cores_number is the number of CPU cores.
# Output: Process and Thread_nums.
# Suppose there are N processes. Process is a list's list and stores N lists. Each list stores command's id.
# Thread_nums is a list storing N elements and each element is the number of threads corresponding to the list in Process.

def GetProcessRuntime(Command_runtime, process_list, thread_num):
    process_runtime = 0
    for command_id in process_list:
        if thread_num == 1:
            process_runtime += Command_runtime[command_id]
        else:
            process_runtime += 1.25 * Command_runtime[command_id] / thread_num

    return process_runtime

def GetMaxProcessId(Command_runtime, Process, Thread_nums):
    max_process_time = 0
    max_process_id = -1
    
    for process_id in range(len(Process)):
        process_runtime = GetProcessRuntime(Command_runtime, Process[process_id], Thread_nums[process_id])
        if runtime > max_time:
            max_process_time = process_runtime
            max_process_id = process_id

    return max_process_id
        
def GetMinProcessId(Command_runtime, Process, Thread_nums):
    min_process_time1 = GetProcessRuntime(Command_runtime, Process[0], Thread_nums[0])
    min_process_id1 = 0
    min_process_time2 = GetProcessRuntime(Command_runtime, Process[1], Thread_nums[1])
    min_process_id2 = 1

    if min_process_time1 > min_process_id2:
        tmp = min_process_time1
        min_process_time1 = min_process_time2
        min_process_time2 = tmp

        tmp = min_process_id1
        min_process_id1 = min_process_id2
        min_process_id2 = tmp

    for process_id in range(2, len(Process)):
        process_runtime = GetProcessRuntime(Command_runtime, Process[process_id], Thread_nums[process_id])
        if process_runtime <= min_process_time1:
            min_process_time2 = min_process_time1
            min_process_id2 = min_process_id1
            min_process_time1 = process_runtime
            min_process_id1 = process_id
        elif process_runtime > min_process_time1 and process_runtime <= min_process_time2:
            min_process_time2 = process_runtime
            min_process_id2 = process_id
    
    return min_process_id1, min_process_id2

def UpdateTotalRuntime(Command_runtime, Process, Thread_nums):
    max_process_time = 0
    for process_id in range(len(Process)):
        process_runtime = GetProcessRuntime(Command_runtime, Process[process_id], Thread_nums[process_id])
        if runtime > max_time:
            max_process_time = process_runtime

    return max_process_time

def LoadBalance(Command_runtime, cores_number):
    Process = []
    Thread_nums = []
    prev_time = 0
    cur_time = 0
    iteration_steps = 0
    # Each process contains one command.
    for i in range(len(Command_runtime)):
        Process_list = []
        Process_list.append(i)
        Process.append(Process_list)
        Thread_nums.append(1)
        prev_time += Command_runtime[i]
    
    # Loop to get the best process-thread running strategy.
    while true:
        if sum(Thread_nums) <= cores_number:
            max_process_id = GetMaxProcessId(Process, Thread_nums)
            Thread_nums[max_process_id] += 1
        else:
            min_process_id1, min_process_id2 = GetMinProcessId(Process, Thread_nums)
            for elem in Process[min_process_id2]:
                Process[min_process_id1].append(elem)
            Thread_nums[min_process_id1] = max(Thread_nums[min_process_id1], Thread_nums[min_process_id2])
            del Process[min_process_id2]
            del Thread_nums[min_process_id2]
        
        prev_time = cur_time
        cur_time = UpdateTotalRuntime(Command_runtime, Process, Thread_nums)
        iteration_steps += 1

        if iteration_steps >= cores_number or abs(cur_time - prev_time) / cur_time < 0.05:
            break

    return Process, Thread_nums
        



