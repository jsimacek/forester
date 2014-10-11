/*
 * Copyright (C) 2010 Jiri Simacek
 *
 * This file is part of forester.
 *
 * forester is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * forester is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with forester.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef NODE_BUILDER_H
#define NODE_BUILDER_H

#include <vector>
#include <cl/code_listener.h>

#include "types.hh"

/**
 * @file nodebuilder.hh
 * NodeBuilder - constructs heap nodes 
 */

/**
 * @brief  Builder of heap nodes
 *
 * Static class wrapping functions creating heap nodes.
 */
struct NodeBuilder {

	/**
	 * @brief  Creates information about heap node
	 *
	 * Creates information about heap node from a Code Listener type. The output
	 * is an array of @e selectors (SelData objects).
	 *
	 * @param[out]  nodeInfo  The output array of SelData objects
	 * @param[in]   type      Properties of the type assigned to the node
	 * @param[in]   name      Name of the selector
	 * @param[in]   offset    Offset of the node in the flat offset space
	 */
	static void buildNode(
		std::vector<SelData>& nodeInfo,
		const cl_type* type,
		int offset = 0,
		const std::string& name = "");


	/**
	 * @brief  Creates information about heap node
	 *
	 * Creates information about heap node from a Code Listener type. The output
	 * is an array of @e flat offsets of fields in the node. Safe for a structure,
	 * the offset is equal to the @p offset parameter; a structure recursively
	 * expands offsets of all of its substructures.
	 *
	 * @param[out]  nodeInfo  The output array of offsets
	 * @param[in]   type      Properties of the type assigned to the node
	 * @param[in]   offset    Offset of the node in the flat offset space
	 */
	static void buildNode(
		std::vector<size_t>& nodeInfo,
		const cl_type* type,
		int offset = 0);


	/**
	 * @brief  Creates information about heap node
	 *
	 * Creates information about heap node from a vector of Code Listener types.
	 * The output is an array of @e flat offsets of fields in the node. Safe for
	 * a structure, the offset is equal to the @p offset parameter; a structure
	 * recursively expands offsets of all of its substructures.
	 *
	 * @param[out]  nodeInfo    The output array of offsets
	 * @param[in]   components  The array of components
	 * @param[in]   offset      The initial offset
	 */
	static void buildNodes(
		std::vector<size_t>& nodeInfo,
		const std::vector<const cl_type*>& components,
		int offset = 0);
};

#endif
