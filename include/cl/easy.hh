/*
 * Copyright (C) 2010 Kamil Dudka <kdudka@redhat.com>
 *
 * This file is part of predator.
 *
 * predator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * predator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with predator.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef H_GUARD_EASY_H
#define H_GUARD_EASY_H

namespace CodeStorage {
    struct Storage;
}

extern void clEasyRun(CodeStorage::Storage &, const char *configString);

/*
 * FIXME: We use the following nonsense in order to keep the static linkage
 *        working.  Moreover a function definition within a header is usually
 *        a really bad idea.  Sooner or later, this should go away!
 */
extern "C" {
    struct plugin_name;
    struct plugin_gcc_version;
    extern int plugin_init(struct plugin_name *, struct plugin_gcc_version *);

#ifndef __CL_IN
    void __cl_easy_stub(void) {
        plugin_init((struct plugin_name *)0, (struct plugin_gcc_version *)0);
    }
#endif
}

#endif /* H_GUARD_EASY_H */
