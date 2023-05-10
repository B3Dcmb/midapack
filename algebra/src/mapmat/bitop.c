/** @file   bitop.c
    @note  Copyright (c) 2010-2012 APC CNRS Université Paris Diderot. This
program is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
more details. You should have received a copy of the GNU General Public License
along with this program; if not, see http://www.gnu.org/licenses/lgpl.html
**/

int is_pow_2(int n) { return ((n & -n) ^ n); }


int pow_2(int k) {
    int n = 1;
    while (k != 0) {
        n = n << 1;
        k--;
    }
    return n;
}


int log_2(int n) {
    int k = 0;
    while (n != 1 && k < 32) {
        n = n >> 1;
        k++;
    }
    return k;
}
