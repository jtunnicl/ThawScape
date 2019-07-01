#ifndef _TIME_FCN_H_
#define _TIME_FCN_H_

class time_fcn {
public:
    time_fcn(int year = 0, int day = 0, int hour = 0, int minute = 0, int end_year = 0, int end_day = 0);
    bool keep_going();
    int get_year();
    int get_day();
    int get_hour();
    int get_minute();
    int get_end_year();
    int get_end_day();
    void increment(int timestep);
    void print();

private:
	int year;
	int day;           // Everything works on 365 Julian Day system for now
	int hour;          // 24 hr clock
	int minute;
	int LocalTime;     // Decimal Hour
	int UT;
	int end_year;      // Model end year
    int end_day;
};

#endif
