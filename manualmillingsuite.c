#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ghcommon.h"

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

int print_menu(int metric_flag)
{
    int choice;

    printf("\n\n\n\t\t*** Main Menu ***\n");
    printf("\n\t<1> Absolute X,Y to Wheel Turns");
    printf("\n\t<2> Wheel turns to Distance");
    printf("\n\t<3> Relative Distance & Angle to Absolute X,Y");
    printf("\n\t<4> Absolute X1,Y1 & X2,Y2 to Distance & Angle");
    printf("\n\t<5> Angle");
    printf("\n\t<6> Radius");
    printf("\n\t<7> Manifold");
    printf("\n\t<8> Decimal Equivalence");
    printf("\n\t<9> Choose Inches or Millimeters");
    if (metric_flag)
        printf(" (Currently Millimeters)");
    else
        printf(" (Currently Inches)");

    printf("\n\n\t<0> Quit");

    do
    {
        choice = get_int("\n\nEnter Selection: ");
    } while (choice < 0 || choice > 9);

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
        printf("\n\n#%d Position\tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", linecount, x, x / 0.03937008, y, y / 0.03937008, wp);
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

void xy_distance_angle(int metric_flag)
{
    double x1,y1,x2,y2,d,a;

    x1 = fabs(get_fraction("\nEnter X1 Position: "));
    y1 = fabs(get_fraction("\nEnter Y1 Position: "));
    x2 = fabs(get_fraction("\nEnter X2 Position: "));
    y2 = fabs(get_fraction("\nEnter Y2 Position: "));

    if (metric_flag)
    {
        x1 = x1 * 0.03937008;
        y1 = y1 * 0.03937008;
        x2 = x2 * 0.03937008;
        y2 = y2 * 0.03937008;
    }

    d = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
    a = rad_to_deg(atan2(y2-y1,x2-x1)) + 90.0;
    if (a < 0.0) a = a + 360.0;

    printf("\n\nDistance:\t%3.4fin %3.2fmm", d, d / 0.03937008);
    printf("\nAngle:   \t%3.1f deg", a);
}

void abs_xy_wheel(double dt, int metric_flag)
{
    double sx, sy;
    char *wp = NULL;

    sx = fabs(get_fraction("\nEnter X Position: "));
    sy = fabs(get_fraction("\nEnter Y Position: "));

    if (metric_flag)
    {
        sx = sx * 0.03937008;
        sy = sy * 0.03937008;
    }

    wp = wheel_position(sx, sy, dt);
    printf("\n\nWheel Position:\tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", sx, sx / 0.03937008, sy, sy / 0.03937008, wp);
    free_malloc(wp);
}

void wheel_distance(double dt, int metric_flag)
{
    double wt, th, d;

    wt = fabs((double)get_int("\nEnter # Wheel Turns: "));
    th = fabs(get_double("\nEnter # of Thousandths: ") / 1000.0);

    wt = wt * dt;
    d = wt + th;

    printf("\n\nDistance:\t%3.4fin %3.2fmm", d, d / 0.03937008);
}

void rel_xy_wheel(double dt, int metric_flag)
{
    double sx, sy, a, ra, d, dx, dy;
    char *wp = NULL;

    sx = fabs(get_fraction("\nEnter Starting X Position: "));
    sy = fabs(get_fraction("\nEnter Starting Y Position: "));

    if (metric_flag)
    {
        sx = sx * 0.03937008;
        sy = sy * 0.03937008;
    }

    wp = wheel_position(sx, sy, dt);
    printf("\n\nStarting at \tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", sx, sx / 0.03937008, sy, sy / 0.03937008, wp);
    free_malloc(wp);

    while (TRUE)
    {
        d = fabs(get_fraction("\n\n\nEnter Distance (0 to Return to Menu): "));
        if (d <= 0.0)
            break;

        if (metric_flag)
        {
            d = d * 0.03937008;
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
        printf("\n\nEnding at \tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", dx, dx / 0.03937008, dy, dy / 0.03937008, wp);
        free_malloc(wp);

        sx = dx;
        sy = dy;
    }
}

void manifold(double dt, int metric_flag)
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

    radius = fabs(get_fraction("\nEnter Radius: "));
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
        radius = radius * 0.03937008;
        centerxoffset = centerxoffset * 0.03937008;
        centeryoffset = centeryoffset * 0.03937008;
    }

    linecount = 0;

    for (a = startangle; a < endangle; a += angleinc)
    {
        radians = PI * a / 180.0;

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
        centerxoffset = centerxoffset * 0.03937008;
        centeryoffset = centeryoffset * 0.03937008;
        inc = inc * 0.03937008;
        tool_radius = tool_radius * 0.03937008;
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

    print_xy_position(x, y, dt, linecount);

    for (step = 0; step < steps; step++)
    {
        x = x + inc * sin(radians);
        y = y - inc * cos(radians);

        linecount++;

        print_xy_position(x, y, dt, linecount);
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
        radius = radius * 0.03937008;
        centerxoffset = centerxoffset * 0.03937008;
        centeryoffset = centeryoffset * 0.03937008;
        arc_length = arc_length * 0.03937008;
    }

    angleinc = arc_length / (PI * radius) * 360.0;

    if (c == 0)
        radius -= tool_radius;
    else
        radius += tool_radius;

    printf("\n\nAngle Increment = %3.1f deg", angleinc);

    linecount = 0;

    for (a = startangle; a <= endangle; a += angleinc)
    {
        radians = deg_to_rad(a);

        x = centerxoffset + radius * sin(radians);
        y = centeryoffset - radius * cos(radians);

        linecount++;

        print_xy_position(x, y, dt, linecount);
    }
}

void print_fraction_double(double v)
{
    fraction fract;
    double i, f;

    f = modf(fabs(v), &i);

    if (v < 0)
        i = -i;

    fract = decimal_to_fraction(f, 1e-6);

    if (fract.n == 0)
    {
        printf("%d", (int)i);
        return;
    }

    if (i <= -1.0 || i >= 1.0)
        printf("%d %d/%d", (int)i, fract.n, fract.d);
    else
        printf("%d/%d", fract.n, fract.d);
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
                printf("  %3.4f : %s (%3.4f)\n", value, sz, value - v);
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

int print_surrounding_fractions(double v, int dnom)
{
    double i, f, fl;
    double om, ol, ex;
    char *st = NULL;

    v = fabs(v);

    if (dnom == 0)
        return FAIL_NUMBER;

    f = modf(v, &i);
    fl = round(f * dnom);

    if (float_compare(fl, 0.0, 1e-6))
        return SUCCESS;

    sprintf_string(&st, "Nearest 1/%d Fractions: ", dnom);

    printf("\n%25s", st);

    free_malloc(st);

    ol = (i + (fl - 1.0) / (double)dnom) - v;
    ex = (i + fl / (double)dnom) - v;
    om = (i + (fl + 1.0) / (double)dnom) - v;

    fraction_int_string(&st, (int)i, (int)(fl - 1.0), dnom);
    print_padded_string(st, 13);
    free_malloc(st);

    sprintf_string(&st, "(%3.4f)", ol);
    print_padded_string(st, 15);
    free_malloc(st);

    fraction_int_string(&st, (int)i, (int)fl, dnom);
    print_padded_string(st, 13);
    free_malloc(st);

    sprintf_string(&st, "(%3.4f)", ex);
    print_padded_string(st, 15);
    free_malloc(st);

    fraction_int_string(&st, (int)i, (int)(fl + 1.0), dnom);
    print_padded_string(st, 13);
    free_malloc(st);

    sprintf_string(&st, "(%3.4f)", om);
    print_padded_string(st, 15);
    free_malloc(st);

    return SUCCESS;
}

void decimal_equivalence(int metric_flag)
{
    double v;
    double f;
    double mm;
    double f64;
    double nf;
    double df;

    int n, i;

    fraction fract;
    char *s = NULL;
    size_t count;

    while (TRUE)
    {
        v = get_fraction("\nEnter a decimal or fractional number (0 to Return to Menu): ");

        if (v == 0.0)
            return;

        if (metric_flag)
            v = v * 0.03937008;

        mm = v / 0.03937008;

        printf("\nInches: %3.4f", v);
        printf("\nMM: %3.2f", mm);
        printf("\nExact Fraction: ");
        print_fraction_double(v);
        printf("\n");

        print_surrounding_fractions(v, 2);
        print_surrounding_fractions(v, 4);
        print_surrounding_fractions(v, 8);
        print_surrounding_fractions(v, 10);
        print_surrounding_fractions(v, 16);
        print_surrounding_fractions(v, 32);
        print_surrounding_fractions(v, 64);
        print_surrounding_fractions(v, 100);
        print_surrounding_fractions(v, 128);

        printf("\n\nStandard Drill Sizes +/- 5%%:\n");
        lookup_drill_size(v);

        printf("\n\n");
    }
}

int main(void)
{
    int choice, metric_flag = 0;
    double dt;

    printf("\n\n\n");
    printf("\t\t\tManual Milling Suite v1.4");
    printf("\n\n");

    print_layout();

    dt = fabs(get_fraction("\nEnter Distance Per Turn of Wheel in Inches: "));

    while (TRUE)
    {
        choice = print_menu(metric_flag);

        printf("\n\n");

        switch (choice)
        {
        case 0:
            return SUCCESS;
        case 1:
            abs_xy_wheel(dt, metric_flag);
            break;
        case 2:
            wheel_distance(dt, metric_flag);
            break;
        case 3:
            rel_xy_wheel(dt, metric_flag);
            break;
        case 4:
            xy_distance_angle(metric_flag);
            break;
        case 5:
            angle(dt, metric_flag);
            break;
        case 6:
            radius(dt, metric_flag);
            break;
        case 7:
            manifold(dt, metric_flag);
            break;
        case 8:
            decimal_equivalence(metric_flag);
            break;
        case 9:
            metric_flag = get_english_or_metric();
            break;
        }
    }
}