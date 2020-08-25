#include <stdio.h>
#include <math.h>

int main()
{
  /* Points (quads) in reference cataloge in counterclockwise direction. */
  double x0=0, x1=2, x2=2, x3=0;
  double y0=0, y1=0, y2=2, y3=2;

  /* Points (quads) in query catalog after matching in CC direction. */
  double ra0=0, ra1=3.46, ra2=1.46, ra3=-2;
  double dec0=0, dec1=2, dec2=5.46, dec3=3.46;

  /* Test Point */
  double x=x1, y=y1;
  double ra=ra1, dec=dec1;


  double X=((x-x0)*ra+(y-y0)*dec)/(ra*ra+dec*dec);
  double Y=(ra*X-(x-x0))/dec;
  double theta=atan(Y/X);
  double s=X/cos(theta);

  printf("X=%g, Y=%g => [theta=%g, s=%g]\n", X, Y, theta*180/3.141592654, s);
}