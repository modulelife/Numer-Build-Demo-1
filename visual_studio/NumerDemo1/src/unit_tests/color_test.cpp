#include "color_test.h"

#include "benchmark.h"
#include <iostream>
#include <string>
#include <numer_visualize.h>
#include <numer_mat.h>
#include <stbi.h>

using namespace numer;






static inline
numer::RGB linearMix(double ratio_1, const numer::RGB& first, const numer::RGB& second)
{
	return numer::RGB{
		(uint8_t)(ratio_1 * first.R + (1.0 - ratio_1) * second.R),
		(uint8_t)(ratio_1 * first.G + (1.0 - ratio_1) * second.G),
		(uint8_t)(ratio_1 * first.B + (1.0 - ratio_1) * second.B)
	};
}

// 1. ����ɽ�ȣ������ɫ���䣬����������Ȼ��
static numer::RGB clrMorningMist(double X) {
    constexpr double t0 = 0.0, t1 = 0.2, t2 = 0.4, t3 = 0.6, t4 = 0.8, t5 = 1.0;

    constexpr numer::RGB
        Clr0{ 0xfb, 0xfc, 0xf9 },  // ���ذ�
        Clr1{ 0xe6, 0xf2, 0xeb },  // ������
        Clr2{ 0xc4, 0xe3, 0xd6 },  // Ϫ����
        Clr3{ 0x94, 0xc7, 0xb0 },  // ާ����
        Clr4{ 0x5a, 0xa5, 0x8e },  // ������
        Clr5{ 0x34, 0x65, 0x54 };  // ɭ����

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else return Clr5;
}

// 2. �����ʸУ��������ɫϵ�����Ƹ߼�֯��
static numer::RGB clrCashmere(double X) {
    constexpr double t0 = 0.0, t1 = 0.15, t2 = 0.35, t3 = 0.55, t4 = 0.75, t5 = 0.9, t6 = 1.0;

    constexpr numer::RGB
        Clr0{ 0xfd, 0xfb, 0xf7 },  // ��Ȼ�ް�
        Clr1{ 0xf8, 0xf2, 0xe7 },  // �װ�
        Clr2{ 0xee, 0xe3, 0xd5 },  // �Һ�
        Clr3{ 0xd7, 0xc9, 0xb9 },  // ������
        Clr4{ 0xba, 0xa8, 0x99 },  // ����ɫ
        Clr5{ 0x9d, 0x8b, 0x7f },  // �̻�
        Clr6{ 0x7a, 0x6c, 0x63 };  // ��Һ�

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else if (X < t6) return linearMix(1.0 - (X - t5) / (t6 - t5), Clr5, Clr6);
    else return Clr6;
}

// 3. ����ˮ�ʣ�����������䣬��ӡ���ɻ���
static numer::RGB clrAqueousBlue(double X) {
    constexpr double t0 = 0.0, t1 = 0.18, t2 = 0.35, t3 = 0.52, t4 = 0.69, t5 = 0.86, t6 = 1.0;

    constexpr numer::RGB
        Clr0{ 0xfd, 0xfe, 0xff },  // ѩ��
        Clr1{ 0xee, 0xf4, 0xfa },  // ���
        Clr2{ 0xd6, 0xe6, 0xef },  // ������
        Clr3{ 0xb3, 0xcc, 0xdc },  // ����
        Clr4{ 0x86, 0xa9, 0xbf },  // ����
        Clr5{ 0x57, 0x82, 0x99 },  // ˮѼ��
        Clr6{ 0x32, 0x55, 0x64 };  // ��ʯ����

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else if (X < t6) return linearMix(1.0 - (X - t5) / (t6 - t5), Clr5, Clr6);
    else return Clr6;
}

// 4. �����ϵ��������ɫ���䣬����֪��
static numer::RGB clrLavenderHarmony(double X) {
    constexpr double t0 = 0.0, t1 = 0.15, t2 = 0.35, t3 = 0.55, t4 = 0.75, t5 = 1.0;

    constexpr numer::RGB
        Clr0{ 0xfc, 0xfb, 0xfd },  // ��Ƥ��
        Clr1{ 0xf5, 0xee, 0xf6 },  // ���ϻ�
        Clr2{ 0xe7, 0xd8, 0xec },  // ޹�²�
        Clr3{ 0xd2, 0xbb, 0xdd },  // ����
        Clr4{ 0xac, 0x8c, 0xbd },  // ÷��ɫ
        Clr5{ 0x79, 0x59, 0x87 };  // ��������

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else return Clr5;
}

// 5. ɰ�ҽ��䣺ů������ɫϵ����ɳĮ����
static numer::RGB clrSandstone(double X) {
    constexpr double t0 = 0.0, t1 = 0.12, t2 = 0.3, t3 = 0.5, t4 = 0.7, t5 = 0.85, t6 = 1.0;

    constexpr numer::RGB
        Clr0{ 0xff, 0xfd, 0xfa },  // ���ع�
        Clr1{ 0xfd, 0xf6, 0xe9 },  // ϸɳ��
        Clr2{ 0xfa, 0xec, 0xd7 },  // ɳ��
        Clr3{ 0xf2, 0xdc, 0xbf },  // ɰ��
        Clr4{ 0xdb, 0xc0, 0xa2 },  // ��ʯ
        Clr5{ 0xb1, 0x98, 0x86 },  // ����
        Clr6{ 0x82, 0x6c, 0x61 };  // ��ľ��

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else if (X < t6) return linearMix(1.0 - (X - t5) / (t6 - t5), Clr5, Clr6);
    else return Clr6;
}

// ������棺��ɫ �� ���� �� ���� �� �� �� ���ٻ� �� ����ɫ
static numer::RGB clrBlueGoldFire(double X) {
    constexpr double t0 = 0.0, t1 = 0.05, t2 = 0.15, t3 = 0.25, t4 = 0.35,
        t5 = 0.45, t6 = 0.55, t7 = 0.65, t8 = 0.75, t9 = 0.85, t10 = 1.0;

    constexpr numer::RGB
        Clr0{ 0x00, 0x00, 0x00 },  // ����
        Clr1{ 0x03, 0x04, 0x17 },  // ��պ�
        Clr2{ 0x08, 0x0c, 0x3a },  // ���
        Clr3{ 0x12, 0x1b, 0x6d },  // Ⱥ����
        Clr4{ 0x1e, 0x42, 0xa3 },  // �ʼ���
        Clr5{ 0x28, 0x77, 0xb8 },  // ������
        Clr6{ 0x31, 0xb0, 0xd1 },  // ����ʯ
        Clr7{ 0xff, 0xd7, 0x00 },  // ����
        Clr8{ 0xff, 0xae, 0x00 },  // ���Ƚ�
        Clr9{ 0xff, 0x8c, 0x00 },  // ������
        Clr10{ 0xff, 0xf9, 0xe6 }; // ů��ɫ

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else if (X < t6) return linearMix(1.0 - (X - t5) / (t6 - t5), Clr5, Clr6);
    else if (X < t7) return linearMix(1.0 - (X - t6) / (t7 - t6), Clr6, Clr7);
    else if (X < t8) return linearMix(1.0 - (X - t7) / (t8 - t7), Clr7, Clr8);
    else if (X < t9) return linearMix(1.0 - (X - t8) / (t9 - t8), Clr8, Clr9);
    else if (X < t10) return linearMix(1.0 - (X - t9) / (t10 - t9), Clr9, Clr10);
    else return Clr10;
}

// �ǿ�������������������ƵĹ���
static numer::RGB clrStellarGold(double X) {
    constexpr double t0 = 0.0, t1 = 0.08, t2 = 0.16, t3 = 0.24, t4 = 0.32,
        t5 = 0.40, t6 = 0.48, t7 = 0.56, t8 = 0.64,
        t9 = 0.72, t10 = 0.80, t11 = 0.88, t12 = 1.0;

    constexpr numer::RGB
        Clr0{ 0x00, 0x00, 0x00 },  // ����
        Clr1{ 0x04, 0x05, 0x20 },  // �����
        Clr2{ 0x0b, 0x10, 0x45 },  // �����
        Clr3{ 0x14, 0x20, 0x70 },  // ��ҹ��
        Clr4{ 0x1f, 0x3d, 0x9f },  // ��ʯ��
        Clr5{ 0x29, 0x6a, 0xc8 },  // ����
        Clr6{ 0x34, 0xa8, 0xdf },  // ������
        Clr7{ 0x42, 0xd4, 0xe8 },  // ����ɫ
        Clr8{ 0x80, 0xe8, 0xd6 },  // ����ʯ
        Clr9{ 0xff, 0xdd, 0x35 },  // ����
        Clr10{ 0xff, 0xb4, 0x20 }, // ���Ƚ�
        Clr11{ 0xff, 0x8c, 0x1a }, // ������
        Clr12{ 0xff, 0xf4, 0xdd }; // ů��ɫ

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else if (X < t6) return linearMix(1.0 - (X - t5) / (t6 - t5), Clr5, Clr6);
    else if (X < t7) return linearMix(1.0 - (X - t6) / (t7 - t6), Clr6, Clr7);
    else if (X < t8) return linearMix(1.0 - (X - t7) / (t8 - t7), Clr7, Clr8);
    else if (X < t9) return linearMix(1.0 - (X - t8) / (t9 - t8), Clr8, Clr9);
    else if (X < t10) return linearMix(1.0 - (X - t9) / (t10 - t9), Clr9, Clr10);
    else if (X < t11) return linearMix(1.0 - (X - t10) / (t11 - t10), Clr10, Clr11);
    else if (X < t12) return linearMix(1.0 - (X - t11) / (t12 - t11), Clr11, Clr12);
    else return Clr12;
}

void ColorTest::run()
{
	stbi::ImageWriter<stbi::format::PNG> writer;
	stbi::ImageLoader loader;


	std::cout << "\n>..ColorTest: begin" << std::endl;
	BENCHMARK_BEGIN(color_test);

	if (loader.load("./image/color/original/gradient.png", stbi::LOAD_GRY))
	{
		mat<uint8_t> img(loader.height(), loader.width(), 0);
		loader.putInto(img.begin());

		writer.writeInto("./image/color/grayscale", img.cbegin(), img.ncols(), img.nrows());

		auto img_num = mat<double>::creat_par(img, GrayScale(0.0, 1.0));

		auto img_re = mat<RGB>::creat_par(img_num, clrStellarGold);

		writer.writeInto("./image/color/newcolor", img_re[0], img_re.ncols(), img_re.nrows());
	}
	
	std::cout << "\n>..ColorTest: end" << std::endl;
	BENCHMARK_END(color_test);
}
