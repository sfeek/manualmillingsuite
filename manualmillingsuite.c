#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ghcommon.h"

#define MM 0.03937008
#define M 3.2808399

double double_mod(double a, double b)
{
    double mod;

    if (a < 0)
        mod = -a;
    else
        mod = a;

    if (b < 0)
        b = -b;

    while (mod >= b)
        mod = mod - b;

    if (a < 0)
        return -mod;

    return mod;
}

int print_main_menu(int metric_flag)
{
    int choice;

    printf("\n\n\n\n\t\t*** Main Menu ***\n");

    printf("\n\t<1>  Distance to Wheel Turns");
    printf("\n\t<2>  Wheel turns to Distance");
    printf("\n\t<3>  Absolute X,Y to Wheel Turns");
    printf("\n\t<4>  Relative Distance & Angle to Absolute X,Y");
    printf("\n\t<5>  Absolute X1,Y1 & X2,Y2 to Distance & Angle");
    printf("\n\t<6>  Angle");
    printf("\n\t<7>  Radius");
    printf("\n\t<8>  Bolt Hole Circle");
    printf("\n\t<9>  Ellipse");
    printf("\n\t<10> Decimal Equivalent");
    printf("\n\t<11> Tap Drill and Clearance Hole Size");
    printf("\n\t<12> Countersink Z Distance");
    printf("\n\t<13> Speeds & Feeds");
    printf("\n\t<14> Choose Inches or Millimeters");
    if (metric_flag)
        printf(" (Currently Millimeters)");
    else
        printf(" (Currently Inches)");

    printf("\n\n\t<0> Quit");

    do
    {
        choice = get_int("\n\nEnter Selection: ");
    } while (choice < 0 || choice > 14);

    return choice;
}

int print_speeds_feeds_menu(int metric_flag)
{
    int choice;

    printf("\n\n\n\n\t\t*** Speed and Feeds Menu ***\n");

    printf("\n\t<1>  Speed");
    printf("\n\t<2>  Feed");
    if (metric_flag)
        printf("\n\t<3>  Surface Meters per Minute");
    else
        printf("\n\t<3>  Surface Feet per Minute");
    if (metric_flag)
        printf("\n\t<4>  Millimeters per Tooth");
    else
        printf("\n\t<4>  Inches per Tooth");
    printf("\n\t<5>  Material Removal Rate");

    printf("\n\n\t<0> Main Menu");

    do
    {
        choice = get_int("\n\nEnter Selection: ");
    } while (choice < 0 || choice > 5);

    return choice;
}

char *wheel_position(double x, double y, double dt)
{
    int xt, yt;
    int xr, yr;
    char *st = NULL;

    xt = (int)trunc(x / dt);
    xr = abs((int)round((x / dt - xt) * dt * 1000));
    yt = (int)trunc(y / dt);
    yr = abs((int)round((y / dt - yt) * dt * 1000));

    sprintf_string(&st, "(%dT + %d, %dT + %d)", xt, xr, yt, yr);

    return st;
}

void print_xy_position(double x, double y, double dt, int linecount)
{
    char *wp = NULL;

    if (x < 0.0 || y < 0.0)
    {
        printf("\n\n#%d Position\tError... Outside of 0,0 bounds!", linecount);
    }
    else
    {
        wp = wheel_position(x, y, dt);
        printf("\n\n#%d Position\tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", linecount, x, x / MM, y, y / MM, wp);
        free_malloc(wp);
    }
}

void print_layout(void)
{
    printf("\n\n\n0,0");
    printf("\n  ---------------");
    printf("\n  |      0      |");
    printf("\n  |             |");
    printf("\n  | 270      90 |");
    printf("\n  |             |");
    printf("\n  |     180     |");
    printf("\n  ---------------");

    printf("\n\n Always reference back left corner as 0,0");

    printf("\n\n\n");
}

int get_english_or_metric(void)
{
    int count, metric_flag;
    char *s = NULL;

    while (TRUE)
    {
        count = get_string(&s, "\n\nEnter <I>nches or <M>illimeters: ");

        if (count == 0)
        {
            free_malloc(s);
            continue;
        }

        if (s[0] == 'I' || s[0] == 'i')
        {
            free_malloc(s);
            metric_flag = FALSE;
            break;
        }

        if (s[0] == 'M' || s[0] == 'm')
        {
            free_malloc(s);
            metric_flag = TRUE;
            break;
        }

        free_malloc(s);
    }

    return metric_flag;
}

int get_inside_outside(void)
{
    int count, inside_outside;
    char *s = NULL;

    while (TRUE)
    {
        count = get_string(&s, "\n\nEnter <I>nside or <O>utside Radius: ");

        if (count == 0)
        {
            free_malloc(s);
            continue;
        }

        if (s[0] == 'I' || s[0] == 'i')
        {
            free_malloc(s);
            inside_outside = FALSE;
            break;
        }

        if (s[0] == 'O' || s[0] == 'o')
        {
            free_malloc(s);
            inside_outside = TRUE;
            break;
        }

        free_malloc(s);
    }

    return inside_outside;
}

int get_lcr(void)
{
    int count, lcr;
    char *s = NULL;

    while (TRUE)
    {
        count = get_string(&s, "\n\nEnter <L>eft <C>enter or <R>ight Side of Angle Line: ");

        if (count == 0)
        {
            free_malloc(s);
            continue;
        }

        if (s[0] == 'L' || s[0] == 'l')
        {
            free_malloc(s);
            lcr = 0;
            break;
        }

        if (s[0] == 'C' || s[0] == 'c')
        {
            free_malloc(s);
            lcr = 1;
            break;
        }

        if (s[0] == 'R' || s[0] == 'r')
        {
            free_malloc(s);
            lcr = 2;
            break;
        }

        free_malloc(s);
    }

    return lcr;
}

void xy_distance_angle(double dt, int metric_flag)
{
    double x1, y1, x2, y2, d, a, mx, my;
    double xd, yd;
    char *wp = NULL;

    x1 = fabs(get_fraction("\nEnter X1 Position: "));
    y1 = fabs(get_fraction("\nEnter Y1 Position: "));
    x2 = fabs(get_fraction("\nEnter X2 Position: "));
    y2 = fabs(get_fraction("\nEnter Y2 Position: "));

    if (metric_flag)
    {
        x1 = x1 * MM;
        y1 = y1 * MM;
        x2 = x2 * MM;
        y2 = y2 * MM;
    }

    mx = (x2 - x1) / 2 + x1;
    my = (y2 - y1) / 2 + y1;
    d = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    a = rad_to_deg(atan2(y2 - y1, x2 - x1)) + 90.0;
    if (a < 0.0)
        a = a + 360.0;

    xd = fabs(x2 - x1);
    yd = fabs(y2 - y1);

    printf("\n\nDistance:\t%3.4fin %3.2fmm", d, d / MM);

    printf("\n\nX Difference:\t%3.4fin %3.2fmm", xd, xd / MM);
    printf("\nY Difference:\t%3.4fin %3.2fmm", yd, yd / MM);

    printf("\n\nAngle:   \t%3.1f deg", a);

    wp = wheel_position(mx, my, dt);
    printf("\n\nMidpoint: \tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", mx, mx / MM, my, my / MM, wp);
    free_malloc(wp);
}

void abs_xy_wheel(double dt, int metric_flag)
{
    double sx, sy;
    char *wp = NULL;

    sx = fabs(get_fraction("\nEnter X Position: "));
    sy = fabs(get_fraction("\nEnter Y Position: "));

    if (metric_flag)
    {
        sx = sx * MM;
        sy = sy * MM;
    }

    wp = wheel_position(sx, sy, dt);
    printf("\n\nWheel Position:\tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", sx, sx / MM, sy, sy / MM, wp);
    free_malloc(wp);
}

void distance_wheel(double dt, int metric_flag)
{
    double d;
    int xt, xr;

    d = fabs(get_fraction("\nEnter Distance: "));

    if (metric_flag)
        d = d * MM;

    xt = (int)trunc(d / dt);
    xr = abs((int)round((d / dt - xt) * dt * 1000));

    printf("\n\nWheel Position:\t%dT + %d", xt,xr);
}

void wheel_distance(double dt, int metric_flag)
{
    double wt, th, d;

    wt = fabs((double)get_int("\nEnter # Wheel Turns: "));
    th = fabs(get_double("\nEnter # of Thousandths: ") / 1000.0);

    wt = wt * dt;
    d = wt + th;

    printf("\n\nDistance:\t%3.4fin %3.2fmm", d, d / MM);
}

void rel_xy_wheel(double dt, int metric_flag)
{
    double sx, sy, a, ra, d, dx, dy;
    char *wp = NULL;

    sx = fabs(get_fraction("\nEnter Starting X Position: "));
    sy = fabs(get_fraction("\nEnter Starting Y Position: "));

    if (metric_flag)
    {
        sx = sx * MM;
        sy = sy * MM;
    }

    wp = wheel_position(sx, sy, dt);
    printf("\n\nStarting at \tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", sx, sx / MM, sy, sy / MM, wp);
    free_malloc(wp);

    while (TRUE)
    {
        d = fabs(get_fraction("\n\n\nEnter Distance (0 to Return to Menu): "));
        if (d <= 0.0)
            break;

        if (metric_flag)
        {
            d = d * MM;
        }

        a = double_mod(get_fraction("\nEnter Angle: "), 360.0);

        ra = PI * a / 180.0;

        dx = sx + d * sin(ra);
        dy = sy - d * cos(ra);

        if (dx < 0.0 || dy < 0.0)
        {
            printf("\n\nError... Outside of 0,0 bounds!");
            continue;
        }

        wp = wheel_position(dx, dy, dt);
        printf("\n\nEnding at \tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", dx, dx / MM, dy, dy / MM, wp);
        free_malloc(wp);

        sx = dx;
        sy = dy;
    }
}

void bolt_hole_circle(double dt, int metric_flag)
{
    double radius;
    double centerxoffset;
    double centeryoffset;
    double startangle;
    double endangle;
    double a;
    double angleinc;
    double x;
    double y;
    double radians;
    int linecount;
    int nh;

    radius = fabs(get_fraction("\nEnter Circle Diameter: ")) / 2.0;
    nh = get_int("\nEnter Number of Holes: ");

    if (nh < 1)
        return;

    angleinc = 360.0 / (double)nh;

    centerxoffset = fabs(get_fraction("\nEnter Center X: "));
    centeryoffset = fabs(get_fraction("\nEnter Center Y: "));
    startangle = double_mod(get_fraction("\nEnter Start Angle: "), 360.0);
    endangle = startangle + 360.0;

    if (metric_flag)
    {
        radius = radius * MM;
        centerxoffset = centerxoffset * MM;
        centeryoffset = centeryoffset * MM;
    }

    linecount = 0;

    for (a = startangle; a < endangle; a += angleinc)
    {
        radians = deg_to_rad(a);

        x = centerxoffset + radius * sin(radians);
        y = centeryoffset - radius * cos(radians);

        linecount++;

        print_xy_position(x, y, dt, linecount);
    }
}

void angle(double dt, int metric_flag)
{
    double tool_radius;
    double centerxoffset;
    double centeryoffset;
    double angle;
    double tangle;
    int steps;
    int step;
    double leftover;
    int lcr;
    double inc;
    double x, y;
    double radians;
    double tradians;
    int linecount = 0;

    centerxoffset = fabs(get_fraction("\nEnter Start X: "));
    centeryoffset = fabs(get_fraction("\nEnter Start Y: "));
    angle = double_mod(get_fraction("\nEnter Angle: "), 360.0);

    inc = fabs(get_fraction("\nEnter Step Distance: "));
    steps = abs(get_int("\nEnter Number of Steps: "));

    tool_radius = fabs(get_fraction("\nEnter Tool Diameter: ") / 2.0);

    lcr = get_lcr();

    if (metric_flag)
    {
        centerxoffset = centerxoffset * MM;
        centeryoffset = centeryoffset * MM;
        inc = inc * MM;
        tool_radius = tool_radius * MM;
    }

    switch (lcr)
    {
    case 0:
        tangle = double_mod(angle - 90, 360.0);
        tradians = deg_to_rad(tangle);
        x = centerxoffset + tool_radius * sin(tradians);
        y = centeryoffset - tool_radius * cos(tradians);
        break;
    case 1:
        x = centerxoffset;
        y = centeryoffset;
        break;
    case 2:
        tangle = double_mod(angle + 90, 360.0);
        tradians = deg_to_rad(tangle);
        x = centerxoffset + tool_radius * sin(tradians);
        y = centeryoffset - tool_radius * cos(tradians);
        break;
    }

    radians = deg_to_rad(angle);

    for (step = 0; step < steps; step++)
    {
        print_xy_position(x, y, dt, linecount);

        x = x + inc * sin(radians);
        y = y - inc * cos(radians);

        linecount++;
    }
}

void radius(double dt, int metric_flag)
{
    double radius;
    double tool_radius;
    double centerxoffset;
    double centeryoffset;
    double startangle;
    double endangle;
    double a;
    double angleinc;
    double x;
    double y;
    double radians;
    double arc_length;
    int linecount;
    int nh;
    int c;

    c = get_inside_outside();

    do
    {
        radius = fabs(get_fraction("\nEnter Radius: "));
        tool_radius = fabs(get_fraction("\nEnter Tool Diameter: ") / 2.0);
    } while ((radius < tool_radius) && (c == 0));

    arc_length = fabs(get_fraction("\nEnter Step Arc Length: "));

    centerxoffset = fabs(get_fraction("\nEnter Center X: "));
    centeryoffset = fabs(get_fraction("\nEnter Center Y: "));
    startangle = double_mod(get_fraction("\nEnter Start Angle: "), 360.0);
    endangle = double_mod(get_fraction("\nEnter End Angle: "), 360.0);

    if (endangle == 0.0)
        endangle = 360.0;

    if (endangle < startangle)
    {
        endangle = endangle + 360.0;
    }

    if (metric_flag)
    {
        radius = radius * MM;
        centerxoffset = centerxoffset * MM;
        centeryoffset = centeryoffset * MM;
        arc_length = arc_length * MM;
        tool_radius = tool_radius * MM;
    }

    angleinc = (arc_length / (PI * radius)) * 360.0;

    if (c == 0)
        radius -= tool_radius;
    else
        radius += tool_radius;

    printf("\n\nAngle Increment:\t%3.1f deg", angleinc);

    linecount = 0;

    for (a = startangle; a <= endangle; a += angleinc)
    {
        radians = deg_to_rad(a);

        x = centerxoffset + radius * sin(radians);
        y = centeryoffset - radius * cos(radians);

        print_xy_position(x, y, dt, linecount);

        linecount++;
    }
}

void print_fraction_double(double v)
{
    fraction fract;
    double i, f;
    char *st = NULL;

    f = modf(fabs(v), &i);

    if (v < 0)
        i = -i;

    fract = decimal_to_fraction(f, 1e-6);

    if (fract.n == 0)
    {
        sprintf_string(&st, "%d", (int)i);
        print_padded_string(st, 13);
        free_malloc(st);
        return;
    }

    if (i <= -1.0 || i >= 1.0)
    {
        sprintf_string(&st, "%d %d/%d", (int)i, fract.n, fract.d);
        print_padded_string(st, 13);
        free_malloc(st);
    }
    else
    {
        sprintf_string(&st, "%d/%d", fract.n, fract.d);
        print_padded_string(st, 13);
        free_malloc(st);
    }
}

void lookup_drill_size(double v)
{
    char line[20];
    char sz[20] = "";
    char last_sz[20] = "";
    double value;
    int flag = 1;

    FILE *file = fopen("drillsize.csv", "r");

    if (file != NULL)
    {
        while (fgets(line, 20, file) != NULL)
        {
            sscanf(line, "%lf,%s", &value, sz);
            if (value >= (v - v * 0.05) && value <= (v + v * 0.05))
            {
                printf("  %3.4f  %s  (%3.4f)\n", value, sz, value - v);
                flag = 0;
            }
        }
    }

    if (file)
        fclose(file);

    if (flag)
        printf("\n  Drill Size Not Found!");

    return;
}

void tap_drill_size(int metric_flag)
{
    double d, p, hp, td, tpi, hd, maxhole;
    int ss;

    if (metric_flag)
    {
        d = fabs(get_fraction("\nEnter Thread Diameter: "));
        p = fabs(get_fraction("\nEnter Thread Pitch: "));
        hp = fabs(get_fraction("\nEnter Thread Depth in %: "));
        hd = fabs(get_fraction("\nEnter Bolt Head Diameter (0 to Skip Clearance Hole Calculation): "));

        td = d - (hp * p) / 76.98;

        printf("\n\nTap Drill Size:\t%3.2f\n", td);

        if (hd > 0.0)
        {
            maxhole = (hd + d) / 2.0;
            printf("\nMin Clearance Hole:\t\t%3.4f", d);
            printf("\nTight Fit Clearance Hole:\t%3.4f", (maxhole - d) * 0.125 + d);
            printf("\nNormal Fit Clearance Hole:\t%3.4f", (maxhole - d) * 0.25 + d);
            printf("\nLoose Fit Clearance Hole:\t%3.4f", (maxhole - d) * 0.50 + d);
            printf("\nMax Clearance Hole:\t\t%3.4f", maxhole);
        }
    }
    else
    {
        printf("\n    Unified Screw Diameter Chart < 1/4in\n");
        for (ss = 0; ss <= 12; ss++)
        {
            printf("\n  #%2d:  %3.3f", ss, 0.060 + ((double)ss * 0.013));
        }

        d = fabs(get_fraction("\n\n\nEnter Thread Diameter: "));
        tpi = fabs(get_fraction("\nEnter Threads per Inch: "));
        hp = fabs(get_fraction("\nEnter Thread Depth in %: "));
        hd = fabs(get_fraction("\nEnter Bolt Head Diameter (0 to Skip Clearance Hole Calculation): "));

        td = d - hp / (76.98 * tpi);

        printf("\n\nTap Drill Size:\t%3.4f\n", td);

        if (hd > 0.0)
        {
            maxhole = (hd + d) / 2.0;
            printf("\nMin Clearance Hole:\t\t%3.4f", d);
            printf("\nTight Fit Clearance Hole:\t%3.4f", (maxhole - d) * 0.125 + d);
            printf("\nNormal Fit Clearance Hole:\t%3.4f", (maxhole - d) * 0.25 + d);
            printf("\nLoose Fit Clearance Hole:\t%3.4f", (maxhole - d) * 0.50 + d);
            printf("\nMax Clearance Hole:\t\t%3.4f", maxhole);
        }
    }
}

void decimal_equivalent(int metric_flag)
{
    double fract = 64.0;
    double v, i, fl;
    double s, e, f;
    double step = 1 / fract;

    while (TRUE)
    {
        v = fabs(get_fraction("\nEnter a decimal or fractional number (0 to Return to Menu): "));

        if (v == 0.0)
            return;

        if (metric_flag)
            v = v * MM;

        printf("\nInches: %3.4f", v);
        printf("\nMM: %3.2f", v / MM);
        printf("\nExact Fraction: ");
        print_fraction_double(v);
        printf("\n\nNearby Fractions");

        f = modf(v, &i);
        fl = round(f * fract);

        s = (fl - 5) / fract + i;
        e = (fl + 5) / fract + i;

        if (s < 0)
            s = 0;

        for (f = s; f <= e; f += step)
        {
            printf("\n  %3.4f  ", f);
            print_fraction_double(f);
            printf("\t(%3.4f)", f - v);
        }

        printf("\n\nStandard Drill Sizes +/- 5%%:\n");
        lookup_drill_size(v);

        printf("\n\n");
    }
}

void countersink(int metric_flag)
{
    double d, a, t, z;

    d = fabs(get_fraction("\nEnter Countersink Diameter: "));
    a = double_mod(fabs(get_fraction("\nEnter Countersink Angle: ")), 360.0);

    if (metric_flag)
        d = d * MM;

    a = deg_to_rad(a);

    z = (d / 2.0) * (1.0 / tan(a / 2.0));

    printf("\n\nZ Depth:\t%3.4fin  %3.2fmm", z, z / MM);
}

void ellipse(double dt, int metric_flag)
{
    double major_radius;
    double minor_radius;
    double foci;
    double rotation_angle;
    double tool_radius;
    double centerxoffset;
    double centeryoffset;
    double startangle;
    double endangle;
    double a;
    double angleinc;
    double x, x0;
    double y, y0;
    double radians;
    double r_radians;
    double arc_length;
    int linecount;
    int c;

    c = get_inside_outside();

    do
    {
        do
        {
            major_radius = fabs(get_fraction("\nEnter Major Radius: "));
            minor_radius = fabs(get_fraction("\nEnter Minor Radius: "));
        } while (major_radius < minor_radius);

        foci = major_radius - sqrt(major_radius * major_radius - minor_radius * minor_radius);
        if (c == 0)
            printf("\n\nMax Tool Diameter:\t%3.4fin  %3.2fmm\n\n", foci, foci / MM);
        tool_radius = fabs(get_fraction("\nEnter Tool Diameter: ") / 2.0);
    } while ((c == 0) && (foci < tool_radius));

    rotation_angle = double_mod(get_fraction("\nEnter Rotation Angle of Minor Axis: "), 360.0);
    arc_length = fabs(get_fraction("\nEnter Step Arc Length: "));

    centerxoffset = fabs(get_fraction("\nEnter Center X: "));
    centeryoffset = fabs(get_fraction("\nEnter Center Y: "));

    startangle = double_mod(get_fraction("\nEnter Start Angle: "), 360.0);
    endangle = double_mod(get_fraction("\nEnter End Angle: "), 360.0);

    if (endangle == 0.0)
        endangle = 360.0;

    if (endangle < startangle)
    {
        endangle = endangle + 360.0;
    }

    if (metric_flag)
    {
        major_radius = major_radius * MM;
        minor_radius = minor_radius * MM;
        centerxoffset = centerxoffset * MM;
        centeryoffset = centeryoffset * MM;
        arc_length = arc_length * MM;
        tool_radius = tool_radius * MM;
    }

    if (c == 0)
    {
        major_radius -= tool_radius;
        minor_radius -= tool_radius;
    }
    else
    {
        major_radius += tool_radius;
        minor_radius += tool_radius;
    }

    linecount = 0;

    r_radians = deg_to_rad(360.0 - double_mod(rotation_angle, 360.0));

    x = 0.0;
    y = 0.0;

    for (a = startangle; a <= endangle; a += 0.001)
    {
        radians = deg_to_rad(a);

        x0 = centerxoffset + (major_radius * cos(radians) * cos(r_radians) - minor_radius * sin(radians) * sin(r_radians));
        y0 = centeryoffset - (major_radius * cos(radians) * sin(r_radians) + minor_radius * sin(radians) * cos(r_radians));

        if (distance(x, y, x0, y0) >= arc_length)
        {
            print_xy_position(x0, y0, dt, linecount);

            x = x0;
            y = y0;

            linecount++;
        }
    }
}

void speed(int metric_flag)
{
    double speed, sfm, d = 0.0;

    if (metric_flag)
    {
        sfm = fabs(get_fraction("\nEnter Surface Meters Per Minute: ")) * 3.2808399;
        do
        {
            d = fabs(get_fraction("\nEnter Tool Diameter: ")) * MM;
        } while (d == 0.0);

        speed = (sfm * 3.82) / d;
        printf("\n\nSpeed (RPM):\t%3.0f", speed);
    }
    else
    {
        sfm = fabs(get_fraction("\nEnter Surface Feet Per Minute: "));
        do
        {
            d = fabs(get_fraction("\nEnter Tool Diameter: "));
        } while (d == 0.0);

        speed = (sfm * 3.82) / d;
        printf("\n\nSpeed (RPM):\t%3.0f", speed);
    }
}

void feed(int metric_flag)
{
    double feed, z, fpt, rpm;

    if (metric_flag)
    {
        rpm = fabs(get_double("\nEnter RPM: "));
        z = fabs((double)get_int("\nEnter # of Teeth: "));
        fpt = fabs(get_fraction("\nEnter Feed per Tooth: "));

        feed = rpm * z * fpt;
        printf("\n\nFeed (MMPM):\t%3.1f", feed);
    }
    else
    {
        rpm = fabs(get_double("\nEnter RPM: "));
        z = fabs((double)get_int("\nEnter # of Teeth: "));
        fpt = fabs(get_fraction("\nEnter Feed per Tooth: "));

        feed = rpm * z * fpt;
        printf("\n\nFeed (IPM):\t%3.1f", feed);
    }
}

void surface(int metric_flag)
{
    double rpm, d, sfm, smm;

    if (metric_flag)
    {
        rpm = fabs(get_double("\nEnter RPM: "));
        d = fabs(get_fraction("\nEnter Tool Diameter: "));

        smm = ((rpm * d * MM) / 3.82) * 0.3048;
        printf("\n\nSurface Meters per Minute:\t%3.3f", smm);
    }
    else
    {
        rpm = fabs(get_double("\nEnter RPM: "));
        d = fabs(get_fraction("\nEnter Tool Diameter: "));

        sfm = (rpm * d) / 3.82;
        printf("\n\nSurface Feet per Minute:\t%3.1f", sfm);
    }
}

void tooth(int metric_flag)
{
    double rpm = 0.0, ipm, mpm, z = 0.0, ipt, mmpt;

    if (metric_flag)
    {
        do
        {
            rpm = fabs(get_fraction("\nEnter RPM: "));
        } while (rpm == 0.0);

        mpm = fabs(get_fraction("\nEnter Feed Rate in Millimeters per Minute: "));

        do
        {
            z = fabs((double)get_int("\nEnter # of Teeth: "));
        } while (z == 0.0);

        mmpt = (mpm / rpm) / z;

        printf("\n\nMillimeters per Tooth:\t%3.2f", mmpt);
    }
    else
    {
        do
        {
            rpm = fabs(get_fraction("\nEnter RPM: "));
        } while (rpm == 0.0);

        ipm = fabs(get_fraction("\nEnter Feed Rate in Inches per Minute: "));

        do
        {
            z = fabs((double)get_int("\nEnter # of Teeth: "));
        } while (z == 0.0);

        ipt = (ipm / rpm) / z;

        printf("\n\nInches per Tooth:\t%3.4f", ipt);
    }
}

void material_removal(int metric_flag)
{
    double ipm, woc, doc;

    if (metric_flag)
    {
        ipm = fabs(get_fraction("\nEnter Feed Rate in Millimeters per Minute: "));
        woc = fabs(get_fraction("\nEnter Width of Cut: "));
        doc = fabs(get_fraction("\nEnter Depth of Cut: "));

        printf("\n\nMaterial Removal in Cubic Millimeters per Minute:\t%3.1f", ipm * woc * doc);
    }
    else
    {
        ipm = fabs(get_fraction("\nEnter Feed Rate in Inches per Minute: "));
        woc = fabs(get_fraction("\nEnter Width of Cut: "));
        doc = fabs(get_fraction("\nEnter Depth of Cut: "));

        printf("\n\nMaterial Removal in Cubic Inches per Minute:\t%3.3f", ipm * woc * doc);
    }
}

int speeds_feeds(int metric_flag)
{
    int choice;

    while (TRUE)
    {
        choice = print_speeds_feeds_menu(metric_flag);

        printf("\n\n");

        switch (choice)
        {
        case 0:
            return SUCCESS;
        case 1:
            speed(metric_flag);
            break;
        case 2:
            feed(metric_flag);
            break;
        case 3:
            surface(metric_flag);
            break;
        case 4:
            tooth(metric_flag);
            break;
        case 5:
            material_removal(metric_flag);
            break;
        }
    }
}

int main(void)
{
    int choice, metric_flag = 0;
    double dt;

    printf("\n\n\n");
    printf("\t\t\tManual Milling Suite v1.85");
    printf("\n\n");

    print_layout();

    dt = fabs(get_fraction("\nEnter Distance Per Turn of Wheel in Inches: "));

    while (TRUE)
    {
        choice = print_main_menu(metric_flag);

        printf("\n\n");

        switch (choice)
        {
        case 0:
            return SUCCESS;
        case 1:
            distance_wheel(dt, metric_flag);
            break;
        case 2:
            wheel_distance(dt, metric_flag);
            break;
        case 3:
            abs_xy_wheel(dt, metric_flag);
            break;
        case 4:
            rel_xy_wheel(dt, metric_flag);
            break;
        case 5:
            xy_distance_angle(dt, metric_flag);
            break;
        case 6:
            angle(dt, metric_flag);
            break;
        case 7:
            radius(dt, metric_flag);
            break;
        case 8:
            bolt_hole_circle(dt, metric_flag);
            break;
        case 9:
            ellipse(dt, metric_flag);
            break;
        case 10:
            decimal_equivalent(metric_flag);
            break;
        case 11:
            tap_drill_size(metric_flag);
            break;
        case 12:
            countersink(metric_flag);
            break;
        case 13:
            speeds_feeds(metric_flag);
            break;
        case 14:
            metric_flag = get_english_or_metric();
            break;
        }
    }
}