/*
 (c) 2011-2015, Vladimir Agafonkin
 SunCalc is a JavaScript library for calculating sun/moon position and light phases.
 https://github.com/mourner/suncalc

 Translated to native C++ by Nikolas Zorko

 Translation notes:
	- I didn't want to use anything outside of the standard library (aka boost) so all dates are currently treated as UnixTime
	- addTime was removed
	- enums are used in place of strings for the TOD data
	- Getting the moon rise/set time is also gonna be a no-go due to treating dates as UnixTime
*/

// sun calculations are based on http://aa.quae.nl/en/reken/zonpositie.html formulas

#pragma once

#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define rad M_PI / 180
#define e_obliquity rad * 23.4397

// date/time constants and conversions

#define dayMS 1000.0 * 60.0 * 60.0 * 24.0
#define J1970 2440588
#define J2000 2451545
#define J0 0.0009

#define SizeOfTimes 6
#define DoubleSizeOfTimes SizeOfTimes * 2

struct SunCoordsData
{
public:
	// declination
	double Declination = 0;

	// right ascension 
	double RightAscension = 0;

	SunCoordsData(double inDeclination, double inRightAscension)
	{
		Declination = inDeclination;
		RightAscension = inRightAscension;
	}
};

struct SunPositionData
{
public:

	double Azimuth;
	double Altitude;

	SunPositionData(double inAzimuth, double inAltitude)
	{
		Azimuth = inAzimuth;
		Altitude = inAltitude;
	}
};

struct MoonCoordsData
{
public:
	double Declination;
	double RightAscension;
	double Distance;

	MoonCoordsData(double inRightAscension, double inDeclination, double inDistance)
	{
		RightAscension = inRightAscension;
		Declination = inDeclination;
		Distance = inDistance;
	}
};

struct MoonPositionData
{
public:
	double Azimuth;
	double Altitude;
	double Distance;
	double ParallacticAngle;

	MoonPositionData(double inAzimuth, double inAltitude, double inDistance, double inParallacticAngle)
	{
		Azimuth = inAzimuth;
		Altitude = inAltitude;
		Distance = inDistance;
		ParallacticAngle = inParallacticAngle;
	}
};

struct MoonIlluminationData
{
public:
	double Fraction;
	double Phase;
	double Angle;

	MoonIlluminationData(double inFraction, double inPhase, double inAngle)
	{
		Fraction = inFraction;
		Phase = inPhase;
		Angle = inAngle;
	}
};

struct MoonRiseData
{
public:
	double SetTime;
	double RiseTime;

	MoonRiseData(double inSetTime, double inRiseTime)
	{
		SetTime = inSetTime;
		RiseTime = inRiseTime;
	}
};

enum SunTOD
{
	sunrise,
	sunriseEnd,
	dawn,
	nauticalDawn,
	nightEnd,
	goldenHourEnd,
	sunset,
	sunsetStart,
	dusk,
	nauticalDusk,
	night,
	goldenHour
};

struct SunCalcTimes
{
public:
	double Angle;
	SunTOD Morning;
	SunTOD Evening;

	SunCalcTimes()
	{
		Angle = 0;
		Morning = SunTOD::sunrise;
		Evening = SunTOD::sunset;
	}

	SunCalcTimes(double inAngle, SunTOD inMorning, SunTOD inEvening)
	{
		Angle = inAngle;
		Morning = inMorning;
		Evening = inEvening;
	}

	double GetAngle() const
	{
		return Angle;
	}
};

struct SunCalcTimeData
{
public:
	double TimeOfDay;
	SunTOD TODType;

	SunCalcTimeData()
	{
		TimeOfDay = 0;
		TODType = SunTOD::sunrise;
	}

	SunCalcTimeData(double inTimeOfDay, SunTOD inTODType)
	{
		TimeOfDay = inTimeOfDay;
		TODType = inTODType;
	}
};

// container for the data returned by getTimes

struct SunCalcTimesContainer
{
private:
	
	int count = 0;
	int maxcount = (DoubleSizeOfTimes);

public:
	double SolarNoon;
	double nadir;

	SunCalcTimeData TimeData[DoubleSizeOfTimes];

	SunCalcTimesContainer()
	{
		SolarNoon = 0;
		nadir = 0;
	}

	SunCalcTimesContainer(double inSolarNoon, double innadir)
	{
		SolarNoon = inSolarNoon;
		nadir = innadir;
	}

	void AddDataToContainer(SunCalcTimeData NewData)
	{
		if (count > maxcount)
			return;

		TimeData[count] = NewData;
		count++;
	}
};

class suncalc
{
public:
	suncalc();
	~suncalc();

private:

	double toJulian(double UnixTime);
	double fromJulian(double JovianTime);
	double toDays(double UnixTime);

	double rightAcension(double l, double b);
	double declination(double l, double b);

	double azimuth(double H, double phi, double dec);
	double altitude(double H, double phi, double dec);

	double siderealTime(double d, double lw);

	double astroRefraction(double h);

	// general sun calculations

	double solarMeanAnomaly(double d);

	double eclipticLongitude(double M);

	SunCoordsData sunCoords(double d);

public:

	SunPositionData getPosition(double UnixTime, double lat, double lon);

private:

	SunCalcTimes times[SizeOfTimes];

	//I would prefer to work with static data, so I've opted to not translate this from JS.
	//void addTime(double Angle, string riseName, string setName);

public:

	// calculations for sun times

	double julianCycle(double d, double lw);

	double approxTransit(double Ht, double lw, double n);
	double solarTransitJ(double ds, double M, double L);

	double hourAngle(double h, double phi, double d);

	// returns set time for the given sun altitude
	double getSetJ(double h, double lw, double phi, double dec, double n, double M, double L);

	// calculates sun times for a given date and latitude/longitude

	SunCalcTimesContainer getTimes(double UnixTime, double Lat, double Lon);

	// moon calculations, based on http://aa.quae.nl/en/reken/hemelpositie.html formulas

	MoonCoordsData moonCoords(double d);

	MoonPositionData getMoonPosition(double UnixTime, double lat, double lng);

	// calculations for illumination parameters of the moon,
	// based on http://idlastro.gsfc.nasa.gov/ftp/pro/astro/mphase.pro formulas and
	// Chapter 48 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.

	MoonIlluminationData getMoonIllumination(double UnixTime);

	double hoursLater(double UnixTime, double h);

	// calculations for moon rise/set times are based on http://www.stargazing.net/kepler/moonrise.html article
	//MoonRiseData getMoonTimes(double UnixTime, double lat, double lon, bool bInUTC = false);
};

