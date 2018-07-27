#include "stdafx.h"
#include "PSO.h"

PSO::PSO(int _max_genom_list, int _var_num, std::vector<double> _varMax, std::vector<double> _varMin) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//dataの初期化
	eliteData(_var_num)
{
	//もらった変数をクラス内変数に格納
	max_genom_list = _max_genom_list;
	var_num = _var_num;
	varMax = _varMax;
	varMin = _varMin;

	for (int i = 0; i < max_genom_list; i++)
	{
		for (int j = 0; j < var_num; j++)
		{
			data[i].x[j] = random(varMin[j], varMax[j]);//遺伝子の初期設定
		}
	}
	prev_data = data;
	calcResult(true);

	displayValues();
}

bool PSO::selection()
{
	int max_num = 0;//最も評価の良い個体の番号
	bool ret = false;

	calc(false);

	for (int i = 0; i < max_genom_list; i++)//ルーレット選択用に評価関数の合計と一番評価の良い番号を取得
	{
		if (data[i].result > data[max_num].result)
			max_num = i;
	}

	eliteData = data[max_num];//最も評価の良い個体を保持
	if (prev_data[max_num].functionValue - eliteData.functionValue != 0)//最も評価の良い個体の変化の監視(デバッグ用)
		ret = true;

	setPosition();
	/*prev_data = data;
	for (int i = 0; i < max_genom_list; i++)
	{
		double selector = random(0.0, 1.0);//乱数を生成
		double needle = 0;//ルーレットの針を生成
		int j = 0;
		for (; ; j++)
		{
			needle += (prev_data[j].result / resultSumValue);//ルーレットの針を乱数の値まで進める
			if (needle > selector)
				break;
			if (j == (max_genom_list - 1))
				break;
		}
		data[i] = prev_data[j];
	}*/
	return ret;
}

/*bool PSO::blxAlphaCrossover()
{
	prev_data = data;

	for (int i = 0; i < max_genom_list; i += 2)//2個ずつ交叉
	{
		for (int j = 0; j < var_num; j++)
		{
			double ave = (data[i].x[j] + data[i + 1].x[j]) / 2;
			double length = std::abs((data[i].x[j] - data[i + 1].x[j]));

			data[i].x[j] = random(ave - length * (1 + alpha * 2) / 2, ave + length * (1 + alpha * 2) / 2);
			data[i + 1].x[j] = random(ave - length * (1 + alpha * 2) / 2, ave + length * (1 + alpha * 2) / 2);
		}
	}
	return true;
}*/

/*bool PSO::mutation()
{
	for (int i = 0; i < max_genom_list; i++)
	{
		if (random(0.0, 1.0) <= individualMutationRate)//個体突然変異率の計算
		{
			for (int j = 0; j < var_num; j++)
			{
				data[i].x[j] = random(varMin[j], varMax[j]);
			}
		}
	}
	return true;
}*/

bool PSO::calc(bool enableDisplay)
{
	calcResult(true);
	for (int i = 0; i < max_genom_list; i++)//評価関数が最小の奴と最大のやつを検索
	{
		if (data[i].result < data[minNum].result)
		{
			minNum = i;
		}
		if (data[i].result > data[maxNum].result)
		{
			maxNum = i;
		}
	}
	//評価関数が最もいいやつを保存
	data[minNum] = eliteData;

	calcResult(true);

	if (enableDisplay)
		displayValues();

	return true;
}

void PSO::setPosition(void)
{
	prev_data = data;
	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].v.size(); j++)
		{
			data[i].v[j] = w * data[i].v[j] + c1 * random(0.0, 1.0)*(data[i].x_pbset[j] - data[i].x[j]) + c2 * random(0.0, 1.0)*(eliteData.x[j]-data[i].x[j]);
		}
	}
	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].x.size(); j++)
		{
			data[i].x[j] += data[i].v[j];
		}
	}
}

bool PSO::calcResult(bool enableSort)
{
	int maxNum = 0;
	double seg;
	
	for (int i = 0; i < max_genom_list; i++)
	{
		data[i].functionValue = std::sin(data[i].x[0] + data[i].x[1]) + std::pow((data[i].x[0] - data[i].x[1]), 2.0) - 1.5*data[i].x[0] + 2.5*data[i].x[1] + 1;//与えられた関数
		if (data[i].functionValue < data[i].functionValuePbset)
		{
			data[i].functionValuePbset = data[i].functionValue;
			data[i].x_pbset = data[i].x;

		}

		/*if (data[maxNum].functionValue < data[i].functionValue)//座標の中で最も関数が大きいやつを検索
		{
			maxNum = i;
		}*/
	}
	seg = data[maxNum].functionValue;//評価関数の切片を与えられた関数が最も大きいやつにセット
	resultSumValue = 0;

	for (int i = 0; i < max_genom_list; i++)
	{
		bool flag = true;
		double coefficient = 0.01;//評価関数用の定数
		for (int j = 0; j < var_num; j++)
		{
			if (data[i].x[j] > varMax[j] || data[i].x[j] < varMin[j])//座標が場外にいるやつの処理
				flag = false;
		}
		data[i].result = std::pow((data[i].functionValue - seg), 2.0);//与えられた関数の値から切片で設定した値を引いて2乗する→与えられた関数の値が小さいやつが強くなる
		//data[i].result = 1 / data[i].functionValue;

		if (!flag)//場外に出たやつの処理
			data[i].result *= coefficient;
		resultSumValue += data[i].result;
	}
	if (enableSort)
		std::sort(data.begin(), data.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });

	return true;
}

int PSO::random(int min, int max)
{
	//乱数の設定
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(engine);
}

double PSO::random(int min, double max)
{
	return random((double)min, max);
}

double PSO::random(double min, int max)
{
	return random(min, (double)max);
}

double PSO::random(double min, double max)
{
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(engine);
}

bool PSO::displayValues()
{
	for (int i = 0; i < max_genom_list; i++)
	{
		for (int j = 0; j < var_num; j++)
		{
			printf_s("%10.7lf,", data[i].x[j]);//デバッグ用
		}
		printf_s(" \t f(x,y)=%10.7lf\t Result=%10.7lf\n", data[i].functionValue, data[i].result);
	}
	return true;
}

PSO::~PSO()
{

}