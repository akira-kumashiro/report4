#include "stdafx.h"
#include "PSO.h"

PSO::PSO(int _max_genom_list, int _var_num, std::vector<double> _varMax, std::vector<double> _varMin) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//data�̏�����
	eliteData(_var_num)
{
	//��������ϐ����N���X���ϐ��Ɋi�[
	varMax = _varMax;
	varMin = _varMin;

	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].x.size(); j++)
		{
			data[i].x[j] = random(varMin[j], varMax[j]);//��`�q�̏����ݒ�
		}
	}
	prev_data = data;
	calcResult(true);

	displayValues();
}

bool PSO::selection()
{
	int max_num = 0;//�ł��]���̗ǂ��̂̔ԍ�
	bool ret = false;

	calc(false);

	for (int i = 0; i < data.size(); i++)//���[���b�g�I��p�ɕ]���֐��̍��v�ƈ�ԕ]���̗ǂ��ԍ����擾
	{
		if (data[i].result > data[max_num].result)
			max_num = i;
	}

	eliteData = data[max_num];//�ł��]���̗ǂ��̂�ێ�
	if (prev_data[max_num].functionValue - eliteData.functionValue != 0)//�ł��]���̗ǂ��̂̕ω��̊Ď�(�f�o�b�O�p)
		ret = true;

	setPosition();
	return ret;
}

void PSO::calc(bool enableDisplay)
{
	calcResult(true);
	for (int i = 0; i < data.size(); i++)//�]���֐����ŏ��̓z�ƍő�̂������
	{
		if (data[i].result < data[minNum].result)
			minNum = i;
	}
	//�]���֐����ł��������ۑ�
	data[minNum] = eliteData;

	calcResult(true);

	if (enableDisplay)
		displayValues();
}

void PSO::setPosition(void)
{
	prev_data = data;
	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].v.size(); j++)//�ʒu�E���x�ݒ�
		{
			data[i].v[j] = w * data[i].v[j] + c1 * random(0.0, 1.0)*(data[i].x_pbset[j] - data[i].x[j]) + c2 * random(0.0, 1.0)*(eliteData.x[j] - data[i].x[j]);
			data[i].x[j] += data[i].v[j];
		}
	}
}

void PSO::calcResult(bool enableSort)
{
	int maxNum = 0;
	double seg;

	for (int i = 0; i < data.size(); i++)
	{
		data[i].functionValue = std::sin(data[i].x[0] + data[i].x[1]) + std::pow((data[i].x[0] - data[i].x[1]), 2.0) - 1.5*data[i].x[0] + 2.5*data[i].x[1] + 1;//�^����ꂽ�֐�
		if (data[i].functionValue < data[i].functionValuePbset)
		{
			data[i].functionValuePbset = data[i].functionValue;
			data[i].x_pbset = data[i].x;

		}
	}
	seg = data[maxNum].functionValue;//�]���֐��̐ؕЂ�^����ꂽ�֐����ł��傫����ɃZ�b�g
	resultSumValue = 0;

	for (int i = 0; i < data.size(); i++)
	{
		bool flag = true;
		double coefficient = 0.01;//�]���֐��p�̒萔
		for (int j = 0; j < data[i].x.size(); j++)
		{
			if (data[i].x[j] > varMax[j] || data[i].x[j] < varMin[j])//���W����O�ɂ����̏���
				flag = false;
		}
		data[i].result = std::pow((data[i].functionValue - seg), 2.0);//�^����ꂽ�֐��̒l����ؕЂŐݒ肵���l��������2�悷�遨�^����ꂽ�֐��̒l����������������Ȃ�

		if (!flag)//��O�ɏo����̏���
			data[i].result *= coefficient;
		resultSumValue += data[i].result;
	}
	if (enableSort)
		std::sort(data.begin(), data.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });
}

int PSO::random(int min, int max)
{
	//�����̐ݒ�
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

void PSO::displayValues()
{
	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].x.size(); j++)
		{
			printf_s("%10.7lf,", data[i].x[j]);//�f�o�b�O�p
		}
		printf_s(" \t f(x,y)=%10.7lf\t Result=%10.7lf\n", data[i].functionValue, data[i].result);
	}
}

PSO::~PSO()
{

}