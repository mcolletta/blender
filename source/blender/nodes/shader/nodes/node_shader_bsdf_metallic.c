/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version. 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2005 Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#include "../node_shader_util.h"

/* **************** OUTPUT ******************** */

static bNodeSocketTemplate sh_node_bsdf_metallic_in[] = {
	{	SOCK_FLOAT,  1, N_("Roughness"),    0.2f, 0.0f, 0.0f, 0.0f,  0.0f, 1.0f, PROP_FACTOR},
	{       SOCK_FLOAT,  1, N_("Anisotropy"),   0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 1.0f},
	{       SOCK_FLOAT,  1, N_("Rotation"),     0.0f, 0.0f, 0.0f, 0.0f,  0.0f, 1.0f, PROP_FACTOR},
	{	SOCK_VECTOR, 1, N_("N"),            0.0f, 0.0f, 0.0f, 1.0f,  0.0f, 100.0f},
	{	SOCK_VECTOR, 1, N_("K"),            0.0f, 0.0f, 0.0f, 1.0f,  0.0f, 100.0f},
	{	SOCK_RGBA,   1, N_("Reflectivity"), 0.8f, 0.8f, 0.8f, 1.0f,  0.0f, 1.0f},
	{	SOCK_RGBA,   1, N_("Edge Tint"),    0.8f, 0.8f, 0.8f, 1.0f,  0.0f, 1.0f},
	{	SOCK_VECTOR, 1, N_("Normal"),	    0.0f, 0.0f, 0.0f, 1.0f, -1.0f, 1.0f, PROP_NONE, SOCK_HIDE_VALUE},
	{       SOCK_VECTOR, 1, N_("Tangent"),      0.0f, 0.0f, 0.0f, 1.0f, -1.0f, 1.0f, PROP_NONE, SOCK_HIDE_VALUE},
	{	-1, 0, ""	}
};

static bNodeSocketTemplate sh_node_bsdf_metallic_out[] = {
	{	SOCK_SHADER, 0, N_("BSDF")},
	{	-1, 0, ""	}
};

static void node_shader_init_metallic(bNodeTree *UNUSED(ntree), bNode *node)
{
	node->custom1 = SHD_METALLIC_ARTISTIC;
}

static void node_shader_update_bsdf_metallic(bNodeTree *UNUSED(ntree), bNode *node)
{
	bNodeSocket *sock;
	int model = node->custom1;

	for (sock = node->inputs.first; sock; sock = sock->next) {
		if (STREQ(sock->name, "N") || STREQ(sock->name, "K")) {
			if (model == SHD_METALLIC_PHYSICAL)
				sock->flag &= ~SOCK_UNAVAIL;
			else
				sock->flag |= SOCK_UNAVAIL;
		}
		else if (STREQ(sock->name, "Reflectivity") || STREQ(sock->name, "Edge Tint")) {
			if (model == SHD_METALLIC_ARTISTIC)
				sock->flag &= ~SOCK_UNAVAIL;
			else
				sock->flag |= SOCK_UNAVAIL;
		}
	}
}

/* node type definition */
void register_node_type_sh_bsdf_metallic(void)
{
	static bNodeType ntype;

	sh_node_type_base(&ntype, SH_NODE_BSDF_METALLIC, "Metallic BSDF", NODE_CLASS_SHADER, 0);
	node_type_compatibility(&ntype, NODE_NEW_SHADING);
	node_type_socket_templates(&ntype, sh_node_bsdf_metallic_in, sh_node_bsdf_metallic_out);
	node_type_size_preset(&ntype, NODE_SIZE_MIDDLE);
	node_type_init(&ntype, node_shader_init_metallic);
	node_type_storage(&ntype, "", NULL, NULL);
	node_type_update(&ntype, node_shader_update_bsdf_metallic, NULL);

	nodeRegisterType(&ntype);
}
