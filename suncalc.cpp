#include "suncalc.h"

// sun times configuration (angle, morning name, evening name)

suncalc::suncalc()
{
	times[0] = SunCalcTimes(-0.833f, SunTOD::sunrise, SunTOD::sunset);
	times[1] = SunCalcTimes(-0.3f, SunTOD::sunriseEnd, SunTOD::sunsetStart);
	times[2] = SunCalcTimes(-6.0f, SunTOD::dawn, SunTOD::dusk);
	times[3] = SunCalcTimes(-12.0f, SunTOD::nauticalDawn, SunTOD::nauticalDusk);
	times[4] = SunCalcTimes(-18.0f, SunTOD::nightEnd, SunTOD::night);
	times[5] = SunCalcTimes(6.0f, SunTOD::goldenHourEnd, SunTOD::goldenHour);

	PI_ = 3.14159265358979323846264338327950288;
	rad = PI_ / 180;
	e_obliquity = rad * 23.4397;

	dayMS = 1000.0 * 60.0 * 60.0 * 24.0;

	J1970 = 2440588;
	J2000 = 2451545;
	J0 = 0.0009;
}

suncalc::~suncalc()
{
}

double suncalc::toJulian(double UnixTime)
{
	return (((UnixTime / dayMS) - 0.5) + J1970);
}

double suncalc::fromJulian(double JovianTime)
{
	return (((JovianTime + 0.5) - J1970) * dayMS);
}

double suncalc::toDays(double UnixTime)
{
	return (toJulian(UnixTime) - J2000);
}

double suncalc::rightAcension(double l, double b)
{
	return atan2(sin(l) * cos(e_obliquity) - tan(b) * sin(e_obliquity), cos(l));
}

double suncalc::declination(double l, double b)
{
	return asin(sin(b) * cos(e_obliquity) + cos(b) * sin(e_obliquity) * sin(l));
}

double suncalc::azimuth(double H, double phi, double dec)
{
	return atan2(sin(H), cos(H) * sin(phi) - tan(dec) * cos(phi));
}

double suncalc::altitude(double H, double phi, double dec)
{
	return asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H));
}

double suncalc::siderealTime(double d, double lw)
{
	return rad * (280.16 + 360.9856235 * d) - lw;
}

double suncalc::astroRefraction(double h)
{
	if (h < 0) // the following formula works for positive altitudes only.
		h = 0; // if h = -0.08901179 a div/0 would occur.

	// formula 16.4 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
	// 1.02 / tan(h + 10.26 / (h + 5.10)) h in degrees, result in arc minutes -> converted to rad:
	return 0.0002967 / tan(h + 0.00312536 / (h + 0.08901179));
}

double suncalc::solarMeanAnomaly(double d)
{
	return rad * (357.5291 + 0.98560028 * d);
}

double suncalc::eclipticLongitude(double M)
{
	double C = rad * (1.9148 * sin(M) + 0.02 * sin(2 * M) + 0.0003 * sin(3 * M)); // equation of center
	double P = rad * 102.9372; // perihelion of the Earth

	return M + C + P + PI_;
}

SunCoordsData suncalc::sunCoords(double d)
{
	double M = solarMeanAnomaly(d);
	double L = eclipticLongitude(M);

	return SunCoordsData(
		declination(L, 0), // dec
		rightAcension(L, 0) // ra
	);
}

SunPositionData suncalc::getPosition(double UnixTime, double lat, double lon)
{
	double lw = rad * -lon;
	double phi = rad * lat;
	double d = toDays(UnixTime);

	SunCoordsData c = sunCoords(d);
	double H = siderealTime(d, lw) - c.RightAscension;

	return SunPositionData(
		azimuth(H, phi, c.Declination),
		altitude(H, phi, c.Declination)
	);
}

double suncalc::julianCycle(double d, double lw)
{
	return round(d - J0 - lw / (2 * PI_));
}

double suncalc::approxTransit(double Ht, double lw, double n)
{
	 return J0 + (Ht + lw) / (2 * PI_) + n;
}

double suncalc::solarTransitJ(double ds, double M, double L)
{
	return J2000 + ds + 0.0053 * sin(M) - 0.0069 * sin(2 * L);
}

double suncalc::hourAngle(double h, double phi, double d)
{
	return acos((sin(h) - sin(phi) * sin(d)) / (cos(phi) * cos(d)));
}

double suncalc::getSetJ(double h, double lw, double phi, double dec, double n, double M, double L)
{
	double w = hourAngle(h, phi, dec);
	double a = approxTransit(w, lw, n);

	return solarTransitJ(a, M, L);
}

SunCalcTimesContainer suncalc::getTimes(double UnixTime, double Lat, double Lon)
{
	double lw = rad * -Lon;
	double phi = rad * Lat;

	double d = toDays(UnixTime);
	double n = julianCycle(d, lw);
	double ds = approxTransit(0, lw, n);

	double M = solarMeanAnomaly(ds);
	double L = eclipticLongitude(M);
	double dec = declination(L, 0);

	double Jnoon = solarTransitJ(ds, M, L);

	int i;
	double Jset;
	double Jrise;

	SunCalcTimesContainer ReturnValue = SunCalcTimesContainer(fromJulian(Jnoon), fromJulian(Jnoon - 0.5));

	for (i = 0; i < SizeOfTimes; i++) 
	{
		SunCalcTimes time = times[i];

		Jset = getSetJ(time.GetAngle() * rad, lw, phi, dec, n, M, L);
		Jrise = Jnoon - (Jset - Jnoon);

		// add two entries per chunk of time data
		ReturnValue.AddDataToContainer(SunCalcTimeData(fromJulian(Jrise), time.Morning));
		ReturnValue.AddDataToContainer(SunCalcTimeData(fromJulian(Jset), time.Evening));
	}

	return ReturnValue;
}

MoonCoordsData suncalc::moonCoords(double d)
{ 
	// geocentric ecliptic coordinates of the moon

	double L = rad * (218.316 + 13.176396 * d); // ecliptic longitude
	double M = rad * (134.963 + 13.064993 * d); // mean anomaly
	double F = rad * (93.272 + 13.229350 * d);  // mean distance

	double l = L + rad * 6.289 * sin(M); // longitude
	double b = rad * 5.128 * sin(F);     // latitude
	double dt = 385001 - 20905 * cos(M);  // distance to the moon in km

	return MoonCoordsData(rightAcension(l, b), declination(l, b), dt);
}

MoonPositionData suncalc::getMoonPosition(double UnixTime, double lat, double lng)
{
	double lw = rad * -lng;
	double phi = rad * lat;
	double d = toDays(UnixTime);

	MoonCoordsData c = moonCoords(d);
	double H = siderealTime(d, lw) - c.RightAscension;
	double h = altitude(H, phi, c.Declination);
		// formula 14.1 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
	double pa = atan2(sin(H), tan(phi) * cos(c.Declination) - sin(c.Declination) * cos(H));

	h = h + astroRefraction(h); // altitude correction for refraction

	return MoonPositionData(azimuth(H, phi, c.Declination), h, c.Distance, pa);
}

MoonIlluminationData suncalc::getMoonIllumination(double UnixTime)
{
	double d = toDays(UnixTime);
	SunCoordsData s = sunCoords(d);
	MoonCoordsData m = moonCoords(d);

	double sdist = 149598000; // distance from Earth to Sun in km

	double phi = acos(sin(s.Declination) * sin(m.Declination) + cos(s.Declination) * cos(m.Declination) * cos(s.RightAscension - m.RightAscension));
	double inc = atan2(sdist * sin(phi), m.Distance - sdist * cos(phi));

	double angle = atan2(cos(s.Declination) * sin(s.RightAscension - m.RightAscension), sin(s.Declination) * cos(m.Declination) - cos(s.Declination) * sin(m.Declination) * cos(s.RightAscension - m.RightAscension));

	return MoonIlluminationData(
		((1 + cos(inc)) / 2), // fraction
		(0.5 + 0.5 * inc * (angle < 0 ? -1 : 1) / PI_), // phase
		angle // angle
	);
}

double suncalc::hoursLater(double UnixTime, double h)
{
	return (UnixTime + h * dayMS / 24);
}

// MoonRiseData suncalc::getMoonTimes(double UnixTime, double lat, double lon, bool bInUTC /*= false*/)
// {
// 	return MoonRiseData(0,0);
// }
