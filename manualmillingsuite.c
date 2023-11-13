#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ghcommon.h"

int print_menu (void)
{
    int choice;

    printf("\n\n\n\t\t*** Main Menu ***\n");
    printf("\n\t<1> Absolute X-Y to Wheel Turns");
    printf("\n\t<2> Absolute X-Y to Relative Angle & Distance");
    printf("\n\t<3> Angle");
    printf("\n\t<4> Radius");
    printf("\n\t<5> Drill Manifold Holes");
    printf("\n\n\t<0> Quit");

    do 
    {
        choice = get_int("\n\nEnter Selection :");
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
        count = get_string(&s, "\n\nEnter <E>nglish or <M>etric: ");

        if (count == 0)
        {
            free_malloc(s);
            continue;
        }

        if (s[0] == 'E' || s[0] == 'e')
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

void abs_xy_wheel(double dt, int metric_flag)
{
    double sx, sy;
    char *wp = NULL;

    sx = get_fraction("\nEnter Starting X Position: ");
    sy = get_fraction("\nEnter Starting Y Position: ");

    if (metric_flag)
    {
        sx = sx * 0.03937008;
        sy = sy * 0.03937008;
    }

    wp = wheel_position(sx, sy, dt);
    printf("\n\nWheel Position:\tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", sx, sx / 0.03937008, sy, sy / 0.03937008, wp);
    free_malloc(wp);
}

void rel_xy_wheel(double dt, int metric_flag)
{
    double sx, sy, a, ra, d, dx, dy;
    char *wp = NULL;

    sx = get_fraction("\nEnter Starting X Position: ");
    sy = get_fraction("\nEnter Starting Y Position: ");

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
        d = get_fraction("\n\n\nEnter Distance (0 to Return to Menu): ");
        if (d <= 0.0)
            break;

        if (metric_flag)
        {
            d = d * 0.03937008;
        }

        a = get_fraction("\nEnter Angle: ");

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

void manifold_wheel (double dt, int metric_flag)
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
    char *wp = NULL;

    radius = fabs(get_fraction("\nEnter Radius: "));
    nh = get_int("\nEnter Number of Holes: ");

    if (nh < 1) return;
    
    angleinc = 360.0 / (double)nh;

    centerxoffset = fabs(get_fraction("\nEnter Center X: "));
    centeryoffset = fabs(get_fraction("\nEnter Center Y: "));
    startangle = get_fraction("\nEnter Start Angle: ");
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

        
        if (x < 0.0 || y < 0.0)
        {
            printf("\n\n#%d Position\tError... Outside of 0,0 bounds!", linecount);
            continue;
        }

        wp = wheel_position(x, y, dt);
        printf("\n\n#%d Position\tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", linecount, x, x / 0.03937008, y, y / 0.03937008, wp);
        free_malloc(wp);
    }
}

void radius_wheel (double dt, int metric_flag)
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
    char *wp = NULL;


    do 
    {
        c = get_int("\nEnter <1> Concave or <2> Convex: ");
    } while (c < 0 || c > 2);

    do
    {
        radius = fabs(get_fraction("\nEnter Radius: "));
        tool_radius = fabs(get_fraction("\nEnter Tool Diameter: ") / 2.0);
    } while ((radius < tool_radius) && (c = 1));

    arc_length = fabs(get_fraction("\nEnter Arc Length: "));
 
    centerxoffset = fabs(get_fraction("\nEnter Center X: "));
    centeryoffset = fabs(get_fraction("\nEnter Center Y: "));
    startangle = get_fraction("\nEnter Start Angle: ");
    endangle = get_fraction("\nEnter End Angle: ");

    if (c==1)
        radius -= tool_radius;
    else
        radius += tool_radius;

    if (metric_flag)
    {
        radius = radius * 0.03937008;
        centerxoffset = centerxoffset * 0.03937008;
        centeryoffset = centeryoffset * 0.03937008;
        arc_length = arc_length * 0.03937008;
    }

    angleinc = arc_length / (PI * radius) * 360.0;
    
    printf("\n\nAngle Increment = %3.1f",angleinc);

    linecount = 0;

    for (a = startangle; a <= endangle; a += angleinc)
    {
        radians = deg_to_rad(a);

        x = centerxoffset + radius * sin(radians);
        y = centeryoffset - radius * cos(radians);

        linecount++;

        if (x < 0.0 || y < 0.0)
        {
            printf("\n\n#%d Position\tError... Outside of 0,0 bounds!", linecount);
            continue;
        }

        wp = wheel_position(x, y, dt);
        printf("\n\n#%d Position\tX = %3.4fin %3.2fmm   Y = %3.4fin %3.2fmm   %s", linecount, x, x / 0.03937008, y, y / 0.03937008, wp);
        free_malloc(wp);
    }
}

int main (void)
{
    int choice, metric_flag;
    double dt;

    printf("\n\n\n");
    printf("\t\t\tManual Milling Suite v1.0");
    printf("\n\n");

    print_layout();
    
    dt = get_fraction("\nEnter Distance Per Turn of Wheel in Inches: ");
    if (dt <= 0.0) return SUCCESS;

    metric_flag = get_english_or_metric();
    
    while (TRUE)
    {
        choice = print_menu();

        printf("\n\n");

        switch (choice)
        {
            case 0:
                return SUCCESS;
            case 1:
                abs_xy_wheel(dt,metric_flag);
                break;
            case 2:
                rel_xy_wheel(dt,metric_flag);
                break;
            case 4:
                radius_wheel(dt,metric_flag);
                break;
            case 5:
                manifold_wheel(dt,metric_flag);
                break;
       }
    }
}