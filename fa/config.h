/*
 * Copyright (C) 2012 Jiri Simacek
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

/**
 * @file config.h
 * various compile-time options
 */

#ifndef CONFIG_H
#define CONFIG_H

/**
 * set abstraction height (default is 1)
 */
#define FA_ABS_HEIGHT						1

/**
 * set reference count tracking treshold (default is 2)
 */
#define FA_REF_CNT_TRESHOLD					2

/**
 * allow to track number of selectors leading towards a given cut-point (default is 0)
 */
#define FA_TRACK_SELECTORS					0

/**
 * allow folding of nested structures (default is 1)
 */
#define FA_ALLOW_FOLDING					1

/**
 * normalize before type 3 learning (default is 0)
 */
#define FA_NORMALIZE_BEFORE_TYPE_3_LEARNING			1

/**
 * unfold type 1 boxes before folding (default is 0)
 */
#define FA_TYPE_1_UNFOLD_HEURISTICS				1

/**
 * should normalization protect type 3 rootpoints (default is FA_ALLOW_FOLDING)
 */
#define FA_KEEP_TYPE_3_CANDIDATES				0

/**
 * overapproximate when folding (default is 0)
 */
#define FA_BOX_APPROXIMATION					1

/**
 * should we restart evry time a new box is encountered (default is 1)
 */
#define FA_RESTART_AFTER_BOX_DISCOVERY				1

/**
 * enable fusion when computing abstraction (default is 1)
 */
#define FA_FUSION_ENABLED					1

#endif /* CONFIG_H */
