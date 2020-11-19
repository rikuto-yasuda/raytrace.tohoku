////////////////////////////////////////////////////////////////////////
// test_model.h
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#ifndef RTC_RAYTRACE_TEST_MODEL_H
#define RTC_RAYTRACE_TEST_MODEL_H

namespace rtc { namespace model {

namespace magnet {

	// test_null model -----------------------------------------
	// test_null���f���ł́A����͑��݂��Ȃ��B
	// ��� 0 ��Ԃ��B
	class test_null_magnet : public basic_magnet_model
	{
	protected: vector getField      ( const vector& pos ) const;
	public:    matrix getDerivativeB( const vector& pos ) const;
	};

	// simpte model ---------------------------------------
	// �P���Ȏ��ꃂ�f�����L�q�B
	// ���̃��f���ł́A����͒P���ȑo�Ɏq����ł���B
	class test_simple : public basic_magnet_model
	{
	public:
		test_simple();

		// test_simple_magnet�͒P���ȑo�Ɏq����ł��邩��A
		// �t�b�g�v�����g�͉�͓I�ɓ������Ƃ��ł���B
		// ���̂��߁Abasic_magnet_model::getFootPrint()��
		// �I�[�o�[���C�h���č��������邱�Ƃ��ł���B
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
	// test_null���f���ł́A�v���Y�}�͑��݂��Ȃ��B
	// ��� 0 ��Ԃ��B
	class test_null_plasma : public basic_plasma_model
	{
	protected:
		double getDensity( const vector& point )   const;
	};

	// test_simple model ---------------------------------------
	// �P���ȃv���Y�}�����f���B
	// ���̃��f���ł́A�n�����ӂɋ��͂ȃv���Y�}�̖�������A
	// x�������� �}�QRe ���ꂽ�ʒu�𒆐S�ɁA���͂ȃv���Y�}�̉򂪑��݂���B
	// �����̃v���Y�}���ɂ��Ȃ�̃f�t�H���������������̂ł���ƍl���č����x���Ȃ��B
	class test_simple : public basic_plasma_model
	{
	protected:
		double getDensity( const vector& point )   const;
	};

}}}// namespace rtc; ---------------------------------------------------

#endif//RTC_RAYTRACE_TEST_MODEL_H
