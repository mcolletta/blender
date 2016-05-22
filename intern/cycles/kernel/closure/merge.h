/*
 * Copyright 2011-2016 Blender Foundation
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

CCL_NAMESPACE_BEGIN

/* Merging
 *
 * This helps with branched path tracing, to avoid sampling identical BSDFs
 * separate, and also saves closures space. Particularly for shadow and emission
 * where we limit the number of closures more. */

ccl_device_noinline void shader_merge_last_closure_without_data(ShaderData *sd)
{
	const int last_closure = ccl_fetch(sd, num_closure) - 1;
	ShaderClosure *sc = ccl_fetch_array(sd, closure, last_closure);

	for(int i = 0; i < last_closure; i++) {
		ShaderClosure *merge_sc = ccl_fetch_array(sd, closure, i);

		if(merge_sc->type == sc->type) {
			merge_sc->weight += sc->weight;
			merge_sc->sample_weight = fabsf(average(merge_sc->weight));
			ccl_fetch(sd, num_closure) -= 1;
			return;
		}
	}
}

ccl_device_noinline void shader_merge_last_closure_with_data(ShaderData *sd)
{
	const int last_closure = ccl_fetch(sd, num_closure) - 1;
	ShaderClosure *sc = ccl_fetch_array(sd, closure, last_closure);

	for(int i = 0; i < last_closure; i++) {
		ShaderClosure *merge_sc = ccl_fetch_array(sd, closure, i);

#ifdef __OSL__
		if(sc->prim || merge_sc->prim)
			continue;
#endif

		if(!(sc->type == merge_sc->type && sc->data0 == merge_sc->data0 && sc->data1 == merge_sc->data1 && sc->data2 == merge_sc->data2))
			continue;

		if(CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
			if(!isequal(sc->N, merge_sc->N))
				continue;
			else if(CLOSURE_IS_BSDF_ANISOTROPIC(sc->type) && !isequal(sc->T, merge_sc->T))
				continue;
		}

		sc->weight += merge_sc->weight;
		sc->sample_weight = fabsf(average(sc->weight));

		ccl_fetch(sd, num_closure) -= 1;
		kernel_assert(sd->num_closure >= 0);
		return;
	}
}

CCL_NAMESPACE_END
