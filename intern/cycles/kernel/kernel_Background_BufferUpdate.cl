/*
 * Copyright 2011-2015 Blender Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "kernel_split.h"

ccl_device_inline void kernel_path_trace_setup(KernelGlobals *kg, ccl_global uint *rng_state, int sample, int x, int y, ccl_addr_space RNG *rng, ccl_addr_space Ray *ray)
{
	float filter_u;
	float filter_v;

	int num_samples = kernel_data.integrator.aa_samples;

	path_rng_init(kg, rng_state, sample, num_samples, rng, x, y, &filter_u, &filter_v);

	/* sample camera ray */
	float lens_u = 0.0f, lens_v = 0.0f;

	if(kernel_data.cam.aperturesize > 0.0f)
		path_rng_2D(kg, rng, sample, num_samples, PRNG_LENS_U, &lens_u, &lens_v);

	float time = 0.0f;

#ifdef __CAMERA_MOTION__
	if(kernel_data.cam.shuttertime != -1.0f)
		time = path_rng_1D(kg, rng, sample, num_samples, PRNG_TIME);
#endif

	camera_sample(kg, x, y, filter_u, filter_v, lens_u, lens_v, time, ray);
}

/*
 * Note on kernel_ocl_path_trace_Background_BufferUpdate kernel.
 * This is the fourth kernel in the ray tracing logic, and the third
 * of the path iteration kernels. This kernel takes care of rays that hit
 * the background (sceneintersect kernel), and for the rays of
 * state RAY_UPDATE_BUFFER it updates the ray's accumulated radiance in
 * the output buffer. This kernel also takes care of rays that have been determined
 * to-be-regenerated.
 *
 * We will empty QUEUE_HITBG_BUFF_UPDATE_TOREGEN_RAYS queue in this kernel
 *
 * Typically all rays that are in state RAY_HIT_BACKGROUND, RAY_UPDATE_BUFFER
 * will be eventually set to RAY_TO_REGENERATE state in this kernel. Finally all rays of ray_state
 * RAY_TO_REGENERATE will be regenerated and put in queue QUEUE_ACTIVE_AND_REGENERATED_RAYS.
 *
 * The input and output are as follows,
 *
 * rng_coop ---------------------------------------------|--- kernel_ocl_path_trace_Background_BufferUpdate ---|--- PathRadiance_coop
 * throughput_coop --------------------------------------|                                                     |--- L_transparent_coop
 * per_sample_output_buffers ----------------------------|                                                     |--- per_sample_output_buffers
 * Ray_coop ---------------------------------------------|                                                     |--- ray_state
 * PathState_coop ---------------------------------------|                                                     |--- Queue_data (QUEUE_HITBG_BUFF_UPDATE_TOREGEN_RAYS)
 * L_transparent_coop -----------------------------------|                                                     |--- Queue_data (QUEUE_ACTIVE_AND_REGENERATED_RAYS)
 * ray_state --------------------------------------------|                                                     |--- Queue_index (QUEUE_HITBG_BUFF_UPDATE_TOREGEN_RAYS)
 * Queue_data (QUEUE_HITBG_BUFF_UPDATE_TOREGEN_RAYS) ----|                                                     |--- Queue_index (QUEUE_ACTIVE_AND_REGENERATED_RAYS)
 * Queue_index (QUEUE_ACTIVE_AND_REGENERATED_RAYS) ------|                                                     |--- work_array
 * parallel_samples -------------------------------------|                                                     |--- PathState_coop
 * end_sample -------------------------------------------|                                                     |--- throughput_coop
 * kg (globals + data) ----------------------------------|                                                     |--- rng_coop
 * rng_state --------------------------------------------|                                                     |--- Ray
 * PathRadiance_coop ------------------------------------|                                                     |
 * sw ---------------------------------------------------|                                                     |
 * sh ---------------------------------------------------|                                                     |
 * sx ---------------------------------------------------|                                                     |
 * sy ---------------------------------------------------|                                                     |
 * stride -----------------------------------------------|                                                     |
 * work_array -------------------------------------------|                                                     |--- work_array
 * queuesize --------------------------------------------|                                                     |
 * start_sample -----------------------------------------|                                                     |--- work_pool_wgs
 * work_pool_wgs ----------------------------------------|                                                     |
 * num_samples ------------------------------------------|                                                     |
 *
 * note on shader_data : shader_data argument is neither an input nor an output for this kernel. It is just filled and consumed here itself.
 * Note on Queues :
 * This kernel fetches rays from QUEUE_HITBG_BUFF_UPDATE_TOREGEN_RAYS queue.
 *
 * State of queues when this kernel is called :
 * At entry,
 * QUEUE_ACTIVE_AND_REGENERATED_RAYS will be filled with RAY_ACTIVE rays
 * QUEUE_HITBG_BUFF_UPDATE_TOREGEN_RAYS will be filled with RAY_UPDATE_BUFFER, RAY_HIT_BACKGROUND, RAY_TO_REGENERATE rays
 * At exit,
 * QUEUE_ACTIVE_AND_REGENERATED_RAYS will be filled with RAY_ACTIVE and RAY_REGENERATED rays
 * QUEUE_HITBG_BUFF_UPDATE_TOREGEN_RAYS will be empty
 */
__kernel void kernel_ocl_path_trace_Background_BufferUpdate(
	ccl_global char *globals,
	ccl_constant KernelData *data,
	ccl_global char *shader_data,
	ccl_global float *per_sample_output_buffers,
	ccl_global uint *rng_state,
	ccl_global uint *rng_coop,                   /* Required for buffer Update */
	ccl_global float3 *throughput_coop,          /* Required for background hit processing */
	PathRadiance *PathRadiance_coop,  /* Required for background hit processing and buffer Update */
	ccl_global Ray *Ray_coop,                    /* Required for background hit processing */
	ccl_global PathState *PathState_coop,        /* Required for background hit processing */
	ccl_global float *L_transparent_coop,        /* Required for background hit processing and buffer Update */
	ccl_global char *ray_state,                  /* Stores information on the current state of a ray */
	int sw, int sh, int sx, int sy, int stride,
	int rng_state_offset_x,
	int rng_state_offset_y,
	int rng_state_stride,
	ccl_global unsigned int *work_array,         /* Denotes work of each ray */
	ccl_global int *Queue_data,                  /* Queues memory */
	ccl_global int *Queue_index,                 /* Tracks the number of elements in each queue */
	int queuesize,                               /* Size (capacity) of each queue */
	int end_sample,
	int start_sample,
#ifdef __WORK_STEALING__
	ccl_global unsigned int *work_pool_wgs,
	unsigned int num_samples,
#endif
#ifdef __KERNEL_DEBUG__
	DebugData *debugdata_coop,
#endif
	int parallel_samples                         /* Number of samples to be processed in parallel */
	)
{
	ccl_local unsigned int local_queue_atomics;
	if(get_local_id(0) == 0 && get_local_id(1) == 0) {
		local_queue_atomics = 0;
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	int ray_index = get_global_id(1) * get_global_size(0) + get_global_id(0);
	if(ray_index == 0) {
		/* We will empty this queue in this kernel */
		Queue_index[QUEUE_HITBG_BUFF_UPDATE_TOREGEN_RAYS] = 0;
	}
	char enqueue_flag = 0;
	ray_index = get_ray_index(ray_index, QUEUE_HITBG_BUFF_UPDATE_TOREGEN_RAYS, Queue_data, queuesize, 1);

#ifdef __COMPUTE_DEVICE_GPU__
	/* If we are executing on a GPU device, we exit all threads that are not required.
	 * If we are executing on a CPU device, then we need to keep all threads active
	 * since we have barrier() calls later in the kernel. CPU devices
	 * expect all threads to execute barrier statement.
	 */
	if(ray_index == QUEUE_EMPTY_SLOT)
		return;
#endif

#ifndef __COMPUTE_DEVICE_GPU__
	if(ray_index != QUEUE_EMPTY_SLOT) {
#endif
		/* Load kernel globals structure and ShaderData strucuture */
		KernelGlobals *kg = (KernelGlobals *)globals;
		ShaderData *sd = (ShaderData *)shader_data;

#ifdef __KERNEL_DEBUG__
		DebugData *debug_data = &debugdata_coop[ray_index];
#endif
		ccl_global PathState *state = &PathState_coop[ray_index];
		PathRadiance *L = L = &PathRadiance_coop[ray_index];
		ccl_global Ray *ray = &Ray_coop[ray_index];
		ccl_global float3 *throughput = &throughput_coop[ray_index];
		ccl_global float *L_transparent = &L_transparent_coop[ray_index];
		ccl_global uint *rng = &rng_coop[ray_index];

#ifdef __WORK_STEALING__
		unsigned int my_work;
		ccl_global float *initial_per_sample_output_buffers;
		ccl_global uint *initial_rng;
#endif
		unsigned int sample;
		unsigned int tile_x;
		unsigned int tile_y;
		unsigned int pixel_x;
		unsigned int pixel_y;
		unsigned int my_sample_tile;

#ifdef __WORK_STEALING__
		my_work = work_array[ray_index];
		sample = get_my_sample(my_work, sw, sh, parallel_samples, ray_index) + start_sample;
		get_pixel_tile_position(&pixel_x, &pixel_y, &tile_x, &tile_y, my_work, sw, sh, sx, sy, parallel_samples, ray_index);
		my_sample_tile = 0;
		initial_per_sample_output_buffers = per_sample_output_buffers;
		initial_rng = rng_state;
#else // __WORK_STEALING__
		sample = work_array[ray_index];
		int tile_index = ray_index / parallel_samples;
		/* buffer and rng_state's stride is "stride". Find x and y using ray_index */
		tile_x = tile_index % sw;
		tile_y = tile_index / sw;
		my_sample_tile = ray_index - (tile_index * parallel_samples);
#endif
		rng_state += (rng_state_offset_x + tile_x) + (rng_state_offset_y + tile_y) * rng_state_stride;
		per_sample_output_buffers += (((tile_x + (tile_y * stride)) * parallel_samples) + my_sample_tile) * kernel_data.film.pass_stride;

		if(IS_STATE(ray_state, ray_index, RAY_HIT_BACKGROUND)) {
			/* eval background shader if nothing hit */
			if(kernel_data.background.transparent && (state->flag & PATH_RAY_CAMERA)) {
				*L_transparent = (*L_transparent) + average((*throughput));
#ifdef __PASSES__
			if(!(kernel_data.film.pass_flag & PASS_BACKGROUND))
#endif
				ASSIGN_RAY_STATE(ray_state, ray_index, RAY_UPDATE_BUFFER);
			}

			if(IS_STATE(ray_state, ray_index, RAY_HIT_BACKGROUND))
			{
#ifdef __BACKGROUND__
				/* sample background shader */
				float3 L_background = indirect_background(kg, state, ray, sd);
				path_radiance_accum_background(L, (*throughput), L_background, state->bounce);
#endif
				ASSIGN_RAY_STATE(ray_state, ray_index, RAY_UPDATE_BUFFER);
			}
		}

		if(IS_STATE(ray_state, ray_index, RAY_UPDATE_BUFFER)) {
			float3 L_sum = path_radiance_clamp_and_sum(kg, L);
			kernel_write_light_passes(kg, per_sample_output_buffers, L, sample);
#ifdef __KERNEL_DEBUG__
			kernel_write_debug_passes(kg, per_sample_output_buffers, state, debug_data, sample);
#endif
			float4 L_rad = make_float4(L_sum.x, L_sum.y, L_sum.z, 1.0f - (*L_transparent));

			/* accumulate result in output buffer */
			kernel_write_pass_float4(per_sample_output_buffers, sample, L_rad);
			path_rng_end(kg, rng_state, *rng);

			ASSIGN_RAY_STATE(ray_state, ray_index, RAY_TO_REGENERATE);
		}

		if(IS_STATE(ray_state, ray_index, RAY_TO_REGENERATE)) {
#ifdef __WORK_STEALING__
			/* We have completed current work; So get next work */
			int valid_work = get_next_work(work_pool_wgs, &my_work, sw, sh, num_samples, parallel_samples, ray_index);
			if(!valid_work) {
				/* If work is invalid, this means no more work is available and the thread may exit */
				ASSIGN_RAY_STATE(ray_state, ray_index, RAY_INACTIVE);
			}
#else
			if((sample + parallel_samples) >= end_sample) {
				ASSIGN_RAY_STATE(ray_state, ray_index, RAY_INACTIVE);
			}
#endif
			if(IS_STATE(ray_state, ray_index, RAY_TO_REGENERATE)) {
#ifdef __WORK_STEALING__
				work_array[ray_index] = my_work;
				/* Get the sample associated with the current work */
				sample = get_my_sample(my_work, sw, sh, parallel_samples, ray_index) + start_sample;
				/* Get pixel and tile position associated with current work */
				get_pixel_tile_position(&pixel_x, &pixel_y, &tile_x, &tile_y, my_work, sw, sh, sx, sy, parallel_samples, ray_index);
				my_sample_tile = 0;

				/* Remap rng_state according to the current work */
				rng_state = initial_rng + ((rng_state_offset_x + tile_x) + (rng_state_offset_y + tile_y) * rng_state_stride);
				/* Remap per_sample_output_buffers according to the current work */
				per_sample_output_buffers = initial_per_sample_output_buffers
											+ (((tile_x + (tile_y * stride)) * parallel_samples) + my_sample_tile) * kernel_data.film.pass_stride;
#else
				work_array[ray_index] = sample + parallel_samples;
				sample = work_array[ray_index];

				/* Get ray position from ray index */
				pixel_x = sx + ((ray_index / parallel_samples) % sw);
				pixel_y = sy + ((ray_index / parallel_samples) / sw);
#endif

				/* initialize random numbers and ray */
				kernel_path_trace_setup(kg, rng_state, sample, pixel_x, pixel_y, rng, ray);

				if(ray->t != 0.0f) {
					/* Initialize throughput, L_transparent, Ray, PathState; These rays proceed with path-iteration*/
					*throughput = make_float3(1.0f, 1.0f, 1.0f);
					*L_transparent = 0.0f;
					path_radiance_init(L, kernel_data.film.use_light_pass);
					path_state_init(kg, state, rng, sample, ray);
#ifdef __KERNEL_DEBUG__
					debug_data_init(debug_data);
#endif
					ASSIGN_RAY_STATE(ray_state, ray_index, RAY_REGENERATED);
					enqueue_flag = 1;
				} else {
					/*These rays do not participate in path-iteration */
					float4 L_rad = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
					/* accumulate result in output buffer */
					kernel_write_pass_float4(per_sample_output_buffers, sample, L_rad);
					path_rng_end(kg, rng_state, *rng);

					ASSIGN_RAY_STATE(ray_state, ray_index, RAY_TO_REGENERATE);
				}
			}
		}
#ifndef __COMPUTE_DEVICE_GPU__
	}
#endif

	/* Enqueue RAY_REGENERATED rays into QUEUE_ACTIVE_AND_REGENERATED_RAYS; These rays
	 * will be made active during next SceneIntersectkernel
	 */
	enqueue_ray_index_local(ray_index, QUEUE_ACTIVE_AND_REGENERATED_RAYS, enqueue_flag, queuesize, &local_queue_atomics, Queue_data, Queue_index);
}