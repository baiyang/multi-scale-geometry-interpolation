#ifndef PAIR__H
#define PAIR__H

class Pair
{
public:
	int a;
	int b;
	Pair()
	{
	}

	Pair(int _a, int _b)
	{
		a = _a;
		b = _b;
	}

	/*** 类字典序比较方式 ***/
	bool operator<(const Pair &A)const 
	{
		if( a < A.a || (a == A.a && b < A.b) )
		{
			return true;
		}


		return false;
	}


	bool operator==(const Pair &A) const
	{
		return a == A.a && b == A.b;
	}
};



#endif //endif
