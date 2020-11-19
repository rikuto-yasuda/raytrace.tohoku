////////////////////////////////////////////////////////////////////////
// test_model.h
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#ifndef RTC_RAYTRACE_TEST_MODEL_H
#define RTC_RAYTRACE_TEST_MODEL_H

namespace rtc { namespace model {

namespace magnet {

	// test_null model -----------------------------------------
	// test_nullモデルでは、磁場は存在しない。
	// 常に 0 を返す。
	class test_null_magnet : public basic_magnet_model
	{
	protected: vector getField      ( const vector& pos ) const;
	public:    matrix getDerivativeB( const vector& pos ) const;
	};

	// simpte model ---------------------------------------
	// 単純な磁場モデルを記述。
	// このモデルでは、磁場は単純な双極子磁場である。
	class test_simple : public basic_magnet_model
	{
	public:
		test_simple();

		// test_simple_magnetは単純な双極子磁場であるから、
		// フットプリントは解析的に導くことができる。
		// そのため、basic_magnet_model::getFootPrint()を
		// オーバーライドして高速化することができる。
		vector getFootPrint(
			const vector& sp,
			double
		) const;

		vector getEquatorPrint(
			const vector& sp,
			double
		) const;

	protected:
		vector getField( const vector& pos ) const;
	};

} namespace plasma {

	// test_null model -----------------------------------------
	// test_nullモデルでは、プラズマは存在しない。
	// 常に 0 を返す。
	class test_null_plasma : public basic_plasma_model
	{
	protected:
		double getDensity( const vector& point )   const;
	};

	// test_simple model ---------------------------------------
	// 単純なプラズマ圏モデル。
	// このモデルでは、地球周辺に強力なプラズマの膜があり、
	// x軸方向に ±２Re 離れた位置を中心に、強力なプラズマの塊が存在する。
	// 現実のプラズマ圏にかなりのデフォルメをかけたものであると考えて差し支えない。
	class test_simple : public basic_plasma_model
	{
	protected:
		double getDensity( const vector& point )   const;
	};

}}}// namespace rtc; ---------------------------------------------------

#endif//RTC_RAYTRACE_TEST_MODEL_H
