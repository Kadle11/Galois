/**
 * @file
 * @section License
 *
 * This file is part of Galois.  Galois is a framework to exploit
 * amorphous data-parallelism in irregular programs.
 *
 * Galois is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, version 2.1 of the
 * License.
 *
 * Galois is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Galois.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * @section Copyright
 *
 * Copyright (C) 2015, The University of Texas at Austin. All rights
 * reserved.
 *
 * @section Description
 *
 * @author Rashid Kaleem<rashid.kaleem@gmail.com>
 */
#ifndef GALOISGPU_OCL_CLKERNEL_H_
#define GALOISGPU_OCL_CLKERNEL_H_

namespace Galois {
namespace OpenCL {
/*************************************************************************
 *
 *************************************************************************/
struct CL_Kernel {
   cl_kernel kernel;
   CL_Device * device;
   size_t global, local;
   double total_time;
   cl_event event;
   char name[256];
   /*************************************************************************
    *
    *************************************************************************/
   explicit CL_Kernel(CL_Device * d) :
         device(d), total_time(0), event(0) {
      global = local = 0;
      kernel = 0;
   }/*************************************************************************
    *
    *************************************************************************/
   explicit CL_Kernel() :
         device(0), total_time(0), event(0) {
      global = local = 0;
      kernel = 0;
   }
   /*************************************************************************
    *
    *************************************************************************/
   CL_Kernel(CL_Device * d, const char * filename, const char * kernel_name, bool from_string) :
         device(d), total_time(0), event(0) {
      if(from_string){
         init_string(filename, kernel_name);
      }else{
         init(filename, kernel_name);
      }
   }
   /*************************************************************************
    *
    *************************************************************************/
   ~CL_Kernel() {
      if (kernel) {
         clReleaseKernel(kernel);
      }
   }
   /*************************************************************************
    *
    *************************************************************************/
   void init(CL_Device * d, const char * filename, const char * kernel_name) {
      device = d;
      strcpy(name, kernel_name);
      global = local = 0;   //Galois::OpenCL::OpenCL_Setup::get_default_device();
      kernel = Galois::OpenCL::getCLContext()->load_kernel(filename, kernel_name, device);
   }
   /*************************************************************************
    *
    *************************************************************************/
   void init(const char * filename, const char * kernel_name) {
      strcpy(name, kernel_name);
      global = local = 0;   //Galois::OpenCL::OpenCL_Setup::get_default_device();
      kernel = Galois::OpenCL::getCLContext()->load_kernel(filename, kernel_name, device);
   }
   /*************************************************************************
    *
    *************************************************************************/
   void init_string(const char * src, const char * kernel_name) {
      strcpy(name, kernel_name);
      global = local = 0;
      kernel = Galois::OpenCL::getCLContext()->load_kernel_string(src, kernel_name, device);
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_work_size(size_t num_items) {
      local = Galois::OpenCL::getCLContext()->workgroup_size(kernel, device);
      global = (size_t) (ceil(num_items / ((double) local)) * local);
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_work_size(size_t num_items, size_t local_size) {
      local = local_size;
      global = (size_t) (ceil(num_items / ((double) local)) * local);
   }
   /*************************************************************************
    *
    *************************************************************************/
   size_t get_default_workgroup_size() {
      return Galois::OpenCL::getCLContext()->workgroup_size(kernel, device);
   }
   /*************************************************************************
    *
    *************************************************************************/
   cl_ulong get_shared_memory_usage() {
      cl_ulong mem_usage;
      CHECK_CL_ERROR(clGetKernelWorkGroupInfo(kernel, device->id(), CL_KERNEL_LOCAL_MEM_SIZE, sizeof(cl_ulong), &mem_usage, NULL),
            "Error: Failed to get shared memory usage for kernel.");
      return mem_usage;
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg(unsigned int index, size_t sz, const void * val) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, index, sz, val), "Arg, compact is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg_shmem(unsigned int index, size_t sz) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, index, sz, nullptr), "Arg, shared-mem is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg_max_shmem(unsigned int index) {
      size_t sz = device->get_device_shared_memory() - get_shared_memory_usage();
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, index, sz, nullptr), "Arg, max-shared-mem is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   size_t get_avail_shmem() {
      size_t sz = device->get_device_shared_memory() - get_shared_memory_usage();
      return sz;
   }
   /*************************************************************************
    *
    *************************************************************************/
   template<typename T>
   void set_arg(unsigned int index, const T &val) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, index, sizeof(cl_mem), &val->device_ptr()), "Arg, compact is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   template<typename T>
   void set_arg_list(const T &val) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val->device_ptr()), "Arg-1, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   template<typename T1, typename T2>
   void set_arg_list(const T1 &val1, const T2 &val2) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1->device_ptr()), "Arg-1/2, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2->device_ptr()), "Arg-2/2, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   template<typename T1, typename T2, typename T3>
   void set_arg_list(const T1 &val1, const T2 &val2, const T3 &val3) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1->device_ptr()), "Arg-1/3, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2->device_ptr()), "Arg-2/3, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3->device_ptr()), "Arg-3/3, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   template<typename T1, typename T2, typename T3, typename T4>
   void set_arg_list(const T1 &val1, const T2 &val2, const T3 &val3, const T4 &val4) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1->device_ptr()), "Arg-1/4, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2->device_ptr()), "Arg-2/4, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3->device_ptr()), "Arg-3/4, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 3, sizeof(cl_mem), &val4->device_ptr()), "Arg-4/4, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   template<typename T1, typename T2, typename T3, typename T4, typename T5>
   void set_arg_list(const T1 &val1, const T2 &val2, const T3 &val3, const T4 &val4, const T5 &val5) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1->device_ptr()), "Arg-1/5, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2->device_ptr()), "Arg-2/5, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3->device_ptr()), "Arg-3/5, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 3, sizeof(cl_mem), &val4->device_ptr()), "Arg-4/5, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 4, sizeof(cl_mem), &val5->device_ptr()), "Arg-5/5, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
   void set_arg_list(const T1 &val1, const T2 &val2, const T3 &val3, const T4 &val4, const T5 &val5, const T6 &val6) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1->device_ptr()), "Arg-1/6, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2->device_ptr()), "Arg-2/6, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3->device_ptr()), "Arg-3/6, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 3, sizeof(cl_mem), &val4->device_ptr()), "Arg-4/6, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 4, sizeof(cl_mem), &val5->device_ptr()), "Arg-5/6, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 5, sizeof(cl_mem), &val6->device_ptr()), "Arg-6/6, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
   void set_arg_list(const T1 &val1, const T2 &val2, const T3 &val3, const T4 &val4, const T5 &val5, const T6 &val6, const T7 & val7) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1->device_ptr()), "Arg-1/7, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2->device_ptr()), "Arg-2/7, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3->device_ptr()), "Arg-3/7, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 3, sizeof(cl_mem), &val4->device_ptr()), "Arg-4/7, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 4, sizeof(cl_mem), &val5->device_ptr()), "Arg-5/7, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 5, sizeof(cl_mem), &val6->device_ptr()), "Arg-6/7, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 6, sizeof(cl_mem), &val7->device_ptr()), "Arg-7/7, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
   void set_arg_list(const T1 &val1, const T2 &val2, const T3 &val3, const T4 &val4, const T5 &val5, const T6 &val6, const T7 & val7, const T8 & val8) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1->device_ptr()), "Arg-1/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2->device_ptr()), "Arg-2/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3->device_ptr()), "Arg-3/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 3, sizeof(cl_mem), &val4->device_ptr()), "Arg-4/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 4, sizeof(cl_mem), &val5->device_ptr()), "Arg-5/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 5, sizeof(cl_mem), &val6->device_ptr()), "Arg-6/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 6, sizeof(cl_mem), &val7->device_ptr()), "Arg-7/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 7, sizeof(cl_mem), &val8->device_ptr()), "Arg-8/8, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg_list_raw(const cl_mem val) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val), "Arg-1, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg_list_raw(const cl_mem &val1, const cl_mem &val2) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1), "Arg-1/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2), "Arg-2/8, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg_list_raw(const cl_mem &val1, const cl_mem &val2, const cl_mem &val3) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1), "Arg-1/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2), "Arg-2/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3), "Arg-3/8, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg_list_raw(const cl_mem &val1, const cl_mem &val2, const cl_mem &val3, const cl_mem &val4) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1), "Arg-1/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2), "Arg-2/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3), "Arg-3/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 3, sizeof(cl_mem), &val4), "Arg-4/8, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg_list_raw(const cl_mem &val1, const cl_mem &val2, const cl_mem &val3, const cl_mem &val4, const cl_mem &val5) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1), "Arg-1/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2), "Arg-2/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3), "Arg-3/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 3, sizeof(cl_mem), &val4), "Arg-4/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 4, sizeof(cl_mem), &val5), "Arg-5/8, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg_list_raw(const cl_mem &val1, const cl_mem &val2, const cl_mem &val3, const cl_mem &val4, const cl_mem &val5, const cl_mem &val6) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1), "Arg-1/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2), "Arg-2/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3), "Arg-3/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 3, sizeof(cl_mem), &val4), "Arg-4/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 4, sizeof(cl_mem), &val5), "Arg-5/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 5, sizeof(cl_mem), &val6), "Arg-6/8, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg_list_raw(const cl_mem &val1, const cl_mem &val2, const cl_mem &val3, const cl_mem &val4, const cl_mem &val5, const cl_mem &val6, const cl_mem & val7) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1), "Arg-1/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2), "Arg-2/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3), "Arg-3/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 3, sizeof(cl_mem), &val4), "Arg-4/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 4, sizeof(cl_mem), &val5), "Arg-5/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 5, sizeof(cl_mem), &val6), "Arg-6/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 6, sizeof(cl_mem), &val7), "Arg-7/8, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_arg_list_raw(const cl_mem &val1, const cl_mem &val2, const cl_mem &val3, const cl_mem &val4, const cl_mem &val5, const cl_mem &val6, const cl_mem & val7, const cl_mem & val8) {
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 0, sizeof(cl_mem), &val1), "Arg-1/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 1, sizeof(cl_mem), &val2), "Arg-2/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 2, sizeof(cl_mem), &val3), "Arg-3/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 3, sizeof(cl_mem), &val4), "Arg-4/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 4, sizeof(cl_mem), &val5), "Arg-5/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 5, sizeof(cl_mem), &val6), "Arg-6/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 6, sizeof(cl_mem), &val7), "Arg-7/8, is NOT set!");
      Galois::OpenCL::CHECK_CL_ERROR(clSetKernelArg(kernel, 7, sizeof(cl_mem), &val8), "Arg-8/8, is NOT set!");
   }
   /*************************************************************************
    *
    *************************************************************************/
   double get_kernel_time() {
      return total_time;
   }
   /*************************************************************************
    *
    *************************************************************************/
   void set_max_threads() {
      set_work_size(local * Galois::OpenCL::getCLContext()->get_device_eu());
   }
   /*************************************************************************
    *
    *************************************************************************/
   void operator()() {
//      fprintf(stderr, "Launching kernel [%s] [%u,%u]\n", name, local, global);
      Galois::OpenCL::CHECK_CL_ERROR(clEnqueueNDRangeKernel(device->command_queue(), kernel, 1, nullptr, &global, &local, 0, nullptr, &event), name);
   }
   /*************************************************************************
    *
    *************************************************************************/
   void run_task() {
      Galois::OpenCL::CHECK_CL_ERROR(clEnqueueTask(device->command_queue(), kernel, 0, nullptr, &event), name);
   }
   /*************************************************************************
    *
    *************************************************************************/
   void operator()(cl_event & e) {
      Galois::OpenCL::CHECK_CL_ERROR(clEnqueueNDRangeKernel(device->command_queue(), kernel, 1, NULL, &global, &local, 0, NULL, &e), name);
   }
   /*************************************************************************
    *
    *************************************************************************/
   void wait() {
      Galois::OpenCL::CHECK_CL_ERROR(clWaitForEvents(1, &event), "Error in waiting for kernel.");
   }
   /*************************************************************************
    *
    *************************************************************************/
   void operator()(size_t num_events, cl_event * events) {
      Galois::OpenCL::CHECK_CL_ERROR(clEnqueueNDRangeKernel(device->command_queue(), kernel, 1, NULL, &global, &local, num_events, events, &event), name);
   }
   /*************************************************************************
    *
    *************************************************************************/
   float last_exec_time() {
#ifdef _GOPT_CL_ENABLE_PROFILING_
      cl_ulong start_time, end_time;
      Galois::OpenCL::CHECK_CL_ERROR(clWaitForEvents(1, &event), "Waiting for kernel");
      Galois::OpenCL::CHECK_CL_ERROR(clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start_time, NULL), "Kernel start time failed.");
      Galois::OpenCL::CHECK_CL_ERROR(clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end_time, NULL), "Kernel end time failed.");
      return (end_time - start_time) / (float) 1.0e9f;
#else
      return 0.0f;
#endif
   }
   /*************************************************************************
    *
    *************************************************************************/
};
}   //namespace OpenCL
} // namespace Galois

#endif /* GALOISGPU_OCL_CLKERNEL_H_ */
