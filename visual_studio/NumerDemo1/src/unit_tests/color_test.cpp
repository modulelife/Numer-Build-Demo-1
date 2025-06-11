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

// 1. 晨雾山谷：柔和绿色渐变，带来宁静自然感
static numer::RGB clrMorningMist(double X) {
    constexpr double t0 = 0.0, t1 = 0.2, t2 = 0.4, t3 = 0.6, t4 = 0.8, t5 = 1.0;

    constexpr numer::RGB
        Clr0{ 0xfb, 0xfc, 0xf9 },  // 晨曦白
        Clr1{ 0xe6, 0xf2, 0xeb },  // 薄雾绿
        Clr2{ 0xc4, 0xe3, 0xd6 },  // 溪流青
        Clr3{ 0x94, 0xc7, 0xb0 },  // 蕨类绿
        Clr4{ 0x5a, 0xa5, 0x8e },  // 松针绿
        Clr5{ 0x34, 0x65, 0x54 };  // 森林绿

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else return Clr5;
}

// 2. 羊绒质感：柔和中性色系，类似高级织物
static numer::RGB clrCashmere(double X) {
    constexpr double t0 = 0.0, t1 = 0.15, t2 = 0.35, t3 = 0.55, t4 = 0.75, t5 = 0.9, t6 = 1.0;

    constexpr numer::RGB
        Clr0{ 0xfd, 0xfb, 0xf7 },  // 天然棉白
        Clr1{ 0xf8, 0xf2, 0xe7 },  // 米白
        Clr2{ 0xee, 0xe3, 0xd5 },  // 灰褐
        Clr3{ 0xd7, 0xc9, 0xb9 },  // 羊绒棕
        Clr4{ 0xba, 0xa8, 0x99 },  // 燕麦色
        Clr5{ 0x9d, 0x8b, 0x7f },  // 烟灰
        Clr6{ 0x7a, 0x6c, 0x63 };  // 深灰褐

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else if (X < t6) return linearMix(1.0 - (X - t5) / (t6 - t5), Clr5, Clr6);
    else return Clr6;
}

// 3. 灰蓝水彩：冷调灰蓝渐变，如印象派画作
static numer::RGB clrAqueousBlue(double X) {
    constexpr double t0 = 0.0, t1 = 0.18, t2 = 0.35, t3 = 0.52, t4 = 0.69, t5 = 0.86, t6 = 1.0;

    constexpr numer::RGB
        Clr0{ 0xfd, 0xfe, 0xff },  // 雪白
        Clr1{ 0xee, 0xf4, 0xfa },  // 冷白
        Clr2{ 0xd6, 0xe6, 0xef },  // 薄雾蓝
        Clr3{ 0xb3, 0xcc, 0xdc },  // 天青
        Clr4{ 0x86, 0xa9, 0xbf },  // 灰蓝
        Clr5{ 0x57, 0x82, 0x99 },  // 水鸭蓝
        Clr6{ 0x32, 0x55, 0x64 };  // 深石板蓝

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else if (X < t6) return linearMix(1.0 - (X - t5) / (t6 - t5), Clr5, Clr6);
    else return Clr6;
}

// 4. 淡雅紫调：柔和紫色渐变，优雅知性
static numer::RGB clrLavenderHarmony(double X) {
    constexpr double t0 = 0.0, t1 = 0.15, t2 = 0.35, t3 = 0.55, t4 = 0.75, t5 = 1.0;

    constexpr numer::RGB
        Clr0{ 0xfc, 0xfb, 0xfd },  // 羊皮白
        Clr1{ 0xf5, 0xee, 0xf6 },  // 淡紫灰
        Clr2{ 0xe7, 0xd8, 0xec },  // 薰衣草
        Clr3{ 0xd2, 0xbb, 0xdd },  // 紫藤
        Clr4{ 0xac, 0x8c, 0xbd },  // 梅紫色
        Clr5{ 0x79, 0x59, 0x87 };  // 深紫罗兰

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else return Clr5;
}

// 5. 砂岩渐变：暖调中性色系，如沙漠景观
static numer::RGB clrSandstone(double X) {
    constexpr double t0 = 0.0, t1 = 0.12, t2 = 0.3, t3 = 0.5, t4 = 0.7, t5 = 0.85, t6 = 1.0;

    constexpr numer::RGB
        Clr0{ 0xff, 0xfd, 0xfa },  // 晨曦光
        Clr1{ 0xfd, 0xf6, 0xe9 },  // 细沙黄
        Clr2{ 0xfa, 0xec, 0xd7 },  // 沙粒
        Clr3{ 0xf2, 0xdc, 0xbf },  // 砂岩
        Clr4{ 0xdb, 0xc0, 0xa2 },  // 赭石
        Clr5{ 0xb1, 0x98, 0x86 },  // 陶土
        Clr6{ 0x82, 0x6c, 0x61 };  // 红木棕

    if (X < t0) return Clr0;
    else if (X < t1) return linearMix(1.0 - (X - t0) / (t1 - t0), Clr0, Clr1);
    else if (X < t2) return linearMix(1.0 - (X - t1) / (t2 - t1), Clr1, Clr2);
    else if (X < t3) return linearMix(1.0 - (X - t2) / (t3 - t2), Clr2, Clr3);
    else if (X < t4) return linearMix(1.0 - (X - t3) / (t4 - t3), Clr3, Clr4);
    else if (X < t5) return linearMix(1.0 - (X - t4) / (t5 - t4), Clr4, Clr5);
    else if (X < t6) return linearMix(1.0 - (X - t5) / (t6 - t5), Clr5, Clr6);
    else return Clr6;
}

// 蓝金火焰：黑色 → 深蓝 → 蓝绿 → 金 → 亮橘黄 → 近白色
static numer::RGB clrBlueGoldFire(double X) {
    constexpr double t0 = 0.0, t1 = 0.05, t2 = 0.15, t3 = 0.25, t4 = 0.35,
        t5 = 0.45, t6 = 0.55, t7 = 0.65, t8 = 0.75, t9 = 0.85, t10 = 1.0;

    constexpr numer::RGB
        Clr0{ 0x00, 0x00, 0x00 },  // 纯黑
        Clr1{ 0x03, 0x04, 0x17 },  // 深空黑
        Clr2{ 0x08, 0x0c, 0x3a },  // 深海蓝
        Clr3{ 0x12, 0x1b, 0x6d },  // 群青蓝
        Clr4{ 0x1e, 0x42, 0xa3 },  // 皇家蓝
        Clr5{ 0x28, 0x77, 0xb8 },  // 天青蓝
        Clr6{ 0x31, 0xb0, 0xd1 },  // 绿松石
        Clr7{ 0xff, 0xd7, 0x00 },  // 纯金
        Clr8{ 0xff, 0xae, 0x00 },  // 亮橙金
        Clr9{ 0xff, 0x8c, 0x00 },  // 火焰橘
        Clr10{ 0xff, 0xf9, 0xe6 }; // 暖白色

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

// 星空蓝金：增加银河蓝到金黄的过渡
static numer::RGB clrStellarGold(double X) {
    constexpr double t0 = 0.0, t1 = 0.08, t2 = 0.16, t3 = 0.24, t4 = 0.32,
        t5 = 0.40, t6 = 0.48, t7 = 0.56, t8 = 0.64,
        t9 = 0.72, t10 = 0.80, t11 = 0.88, t12 = 1.0;

    constexpr numer::RGB
        Clr0{ 0x00, 0x00, 0x00 },  // 纯黑
        Clr1{ 0x04, 0x05, 0x20 },  // 宇宙黑
        Clr2{ 0x0b, 0x10, 0x45 },  // 深空蓝
        Clr3{ 0x14, 0x20, 0x70 },  // 午夜蓝
        Clr4{ 0x1f, 0x3d, 0x9f },  // 宝石蓝
        Clr5{ 0x29, 0x6a, 0xc8 },  // 钴蓝
        Clr6{ 0x34, 0xa8, 0xdf },  // 天青蓝
        Clr7{ 0x42, 0xd4, 0xe8 },  // 蓝绿色
        Clr8{ 0x80, 0xe8, 0xd6 },  // 绿松石
        Clr9{ 0xff, 0xdd, 0x35 },  // 纯金
        Clr10{ 0xff, 0xb4, 0x20 }, // 亮橙金
        Clr11{ 0xff, 0x8c, 0x1a }, // 火焰橘
        Clr12{ 0xff, 0xf4, 0xdd }; // 暖白色

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
