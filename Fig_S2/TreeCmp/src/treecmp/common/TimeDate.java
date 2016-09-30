/** This file is part of TreeCmp, a tool for comparing phylogenetic trees
    using the Matching Split distance and other metrics.
    Copyright (C) 2011,  Damian Bogdanowicz

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

package treecmp.common;

import java.text.SimpleDateFormat;
import java.util.Calendar;

public class TimeDate {

    public static final String DATE_FORMAT_NOW="yyyy-MM-dd HH:mm:ss";

    public static String now()
    {

        Calendar cal=Calendar.getInstance();
        SimpleDateFormat sdf=new SimpleDateFormat(DATE_FORMAT_NOW);
        return sdf.format(cal.getTime());
        
    }




}
