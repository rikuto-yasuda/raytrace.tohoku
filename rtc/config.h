//////////////////////////////////////////////////////////////////////
// config.h
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#ifndef RTC_RAYTRACE_CORE_CONFIG_H
#define RTC_RAYTRACE_CORE_CONFIG_H

//
// raytrace�̒����\�ȃp�����[�^�͂����ŋL�q�B
//

//
// ���x���z�A������z���v�Z����Ƃ��̕������[�g���P�ʂŎw�肷��B
//
#define RTC_DERIVATIVE_DISTANCE 1.0


//
// basic_magnet_model::getFootPrint()���A
// �ߋ��̌v�Z���ʂ�ێ�������������}��ꍇ�ɒ�`����B
// ��`���ꂽ���������ߋ��̌v�Z���ʂ�ێ�����B
//
#define RTC_BASIC_MAGNET_MODEL_STORE_PAST 5


//
// cavity �ɂ��v���Y�}���x�ւ̉e����^����ꍇ�A��`����B
// ������`���Ȃ����Ƃɂ���āAbasic_plasma_model::operator ()�̌��ʂ�
// �v���Y�}�L���r�e�B���܂܂Ȃ��A���̒l��Ԃ��B
// ���̌��ʁA�ق�̏��������̍������������ނ��Ƃ��ł���B
// ����`�̏�Ԃł� cavity�N���X�̃C���X�^���X�����邱�Ƃ͂ł��邵
// setCavity()���Ăяo�����Ƃ��\�����A�v���Y�}�ɂ͑S���e�����Ȃ��B
//
#define RTC_ENABLE_PLASMA_CAVITY


//
// raytrace���ɁA�p�����[�^�w��̎��Ԃ�菭�Ȃ��o�ߎ��Ԃ��K�v��
// ���f���ꂽ�����G���[�Ɣ��f���A��O���Ȃ���悤�ɂ���B
// ���f���̕s��ȂǂŖ������[�v�Ɋׂ�ꍇ�A�L���ɂ���B
//
#define RTC_RAYTRACE_ENABLE_EXCEPTION_WHEN_TIMESTEP_UNDERFLOW


// -----------------------------------------------
// �ȉ��A�f�o�b�O���ɂ̂ݗL���ɂ���B
// ��f�o�b�O�����L���ɂ���ꍇ��
// ifndef�߂��R�����g�A�E�g����B
#ifndef NDEBUG

//
// rtc::clearNaN()��L���ɂ���ꍇ�A��`����B
//
#define RTC_RAYTRACE_ENABLE_CLEARNAN


//
// rtc::clearNaN()���L���ł��������A
// ���ۂ�NaN���N���A�����Ƃ��Ƀ��|�[�g�\������ꍇ�ɒ�`����B
//
#define RTC_RAYTRACE_LOGGING_CLEARNAN


//
// ray::calc_dt()�ŁA�p�����[�^�w��̎��Ԕ͈͊O�ɓ��B
// �����Ƃ��Ƀ��O���Ƃ�ꍇ�A��`����B
//
#define RTC_RAYTRACE_RAY_LOGS_OUT_OF_TIME_RANGE


//
// �f�o�b�O�p�̃��O���o�͂���ꍇ�A��`����B
//
#define RTC_RAYTRACE_ENABLE_LOG


//
// ���ڍׂȃf�o�b�O���O���o�͂���ꍇ�A��`����B
//
// #define RTC_RAYTRACE_ENABLE_DETAIL_LOG

#endif//NDEBUG -----------------------------------


#endif//RTC_RAYTRACE_CORE_CONFIG_H
