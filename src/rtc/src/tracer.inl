////////////////////////////////////////////////////////////////////////
// tracer.inl
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
namespace rtc {

// ���C�g���[�V���O�����s���A����i�� //////////////////////////////////
double ray::take_a_step()
{
	// ��������x�N�g�������o���B
	vector& r  = m_rk.first;
	vector& k  = m_rk.second;

#ifndef NDEBUG
	intermediate i;
	update_intermediate(i,r,k);

	assert( i == m_im );
#endif

	const double dwdg = 1.0/calc_dGdw(m_im,r,k);
	const vector dgdr = calc_dGdr(m_im,r,k);
	const vector dgdk = calc_dGdk(m_im,r,k);

	m_drk.first  = -dgdk*dwdg; // -drdt
	m_drk.second =  dgdr*dwdg; //  dkdt

	// calc_dt() �̒��� update_intermediate()���Ăяo����
	// m_im ���X�V�����B
	const double dt = calc_dt(
		m_rk, m_drk,
		m_im
	);

	// ���ʂ��i�[���āA�����Ԃ��B
	m_drk.first  *= dt;
	m_drk.second *= dt;
	r += m_drk.first;
	k += m_drk.second;

	if (r(2)<0)
	{
		r(2)=0;
		k(2)=-k(2);
	}

	// �V�����_�ł̌��̏�Ԃ��`�F�b�N����B
	checkState(m_im,r,k);

	return dt;
}

}// namespace rtc
