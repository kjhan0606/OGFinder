"""Guided science recipes and immutable, verified research checkpoints.

This is a guided workflow catalogue, not an autonomous science validator.
Only explicitly enumerated regular output files are copied; inputs are references.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import tempfile
from datetime import datetime, timezone
import uuid


def node(key, title, detail, action=None, optional=False):
    return dict(id=key, title=title, detail=detail, action=action or [], optional=optional)


AUDIT = node('audit', '입력·교정·오차 확인',
             '단위, 영점, WCS, variance, 마스크, 포화, 처리 이력을 확인합니다. 원시 자료는 기구별 교정이 필요합니다.')
EXPORT = node('results', '결과 검토·세션 저장',
              '측정 표·오차·잔차·검증 결과를 확인하고 산출물을 연결한 뒤 연구 세션을 저장합니다. 검토 완료는 과학 검증 통과와 다릅니다.')


def recipe(key, title, inputs, outputs, steps):
    all_steps = [node('load', '자료 불러오기', inputs, ['load']), dict(AUDIT)] + steps + [dict(EXPORT)]
    all_steps[-1]['action'] = ['checkpoint']
    if not key.endswith('_spectrum'):
        all_steps[1]['action'] = ['ds10core', 'workflow-audit']
        all_steps[1]['web_action'] = ['ds10_workflow', 'audit']
        all_steps[1]['detail'] += ' 실행하면 유효/결함 픽셀, 평균·산포와 헤더를 보고합니다. 산포는 sky 잡음 추정값이 아닙니다.'
    # The guided path is linear; conditions stay visible rather than silently omitted.
    return dict(id=key, title=title, version=1, inputs=inputs, outputs=outputs,
                steps=all_steps, edges=[[a['id'], b['id']] for a, b in zip(all_steps, all_steps[1:])])


def recipes():
    star_steps = [
        node('detect', '별 검출·초기 목록', '혼잡도와 포화 별을 확인하고 초기 별 목록을 만듭니다.', ['extract', 'extract']),
        node('aperture', '초기 구경 측광', 'PSF 별 선택에 필요한 초기 등급과 하늘 추정값을 확인합니다.', ['daophot', 'phot']),
        node('pickpsf', 'PSF 별 선택', '고립된 별과 이웃 영향을 검토합니다.', ['daophot', 'pickpsf']),
        node('psf', 'PSF 모형 구성', '공간 변화와 PSF 잔차를 확인합니다.', ['daophot', 'psf']),
        node('fit', 'PSF 측광·이웃 제거', '혼잡한 별을 적합하고 잔차 재검출을 검토합니다.', ['daophot', 'fit']),
        node('calibrate', '구경 보정·측광 교정', '구경 보정 후 영점·색항·오차를 확인합니다. 영점이 없으면 기기등급으로 표시합니다.', ['daophot', 'apcorr']),
        node('match', '노출·밴드 간 대응', '같은 별을 대응합니다. 단일 영상이면 생략 사유를 기록합니다.', ['daophot', 'match'], True),
        node('complete', '인공별 회수·완전도', '등급과 중심 거리별로 최종 품질 컷까지 회수합니다.', ['daophot', 'arttest']),
        node('cmd', 'CMD·성단 멤버 검토', '두 밴드 이상일 때 CMD를 생성합니다. 검출 별 전체를 멤버로 확정하지 않습니다.', ['daophot', 'cmd'], True),
    ]
    result = [recipe('globular_cluster', '구상성단 · 별 측광',
        '교정된 성단 영상, 오차/마스크, 영점; CMD는 두 밴드 이상.',
        '별 측광 표, 오차·완전도, 잔차 영상, 가능한 경우 CMD', star_steps)]
    result.append(recipe('open_cluster', '산개성단 · 측광과 멤버',
        '교정된 성단 영상과 밴드 정보; 멤버 판정에는 고유운동·시차 등 추가 자료.',
        '측광 표, CMD, 근거를 명시한 멤버 후보 및 완전도', [
        node('detect', '별 검출·혼잡도 판정', '혼잡도에 따라 구경 또는 PSF 측광 경로를 선택합니다.', ['extract', 'extract']),
        node('phot', '측광 경로 설정·실행', '분리된 별은 구경 측광, 혼잡한 영역은 PSF 적합을 선택하고 근거를 기록합니다.', ['daophot', 'run']),
        node('calibrate', '구경 보정·영점 확인', '측광계와 색항, 노출 간 교정을 확인합니다.', ['daophot', 'apcorr']),
        node('match', '밴드·측성 카탈로그 대응', '다중 밴드 또는 고유운동·시차 자료를 대응하고 선택함수를 기록합니다.', ['xmatch', 'match']),
        node('membership', '멤버십·소광 검토', '관측 오차와 전경·배경을 포함한 멤버 모형이 필요합니다. 은하단 red-sequence 기능으로 대체하지 않습니다.'),
        node('complete', '인공별 회수', '등급·위치별 검출 및 측광 편향을 평가합니다.', ['daophot', 'arttest']),
        node('cmd', 'CMD 및 해석 범위', '나이·거리 적합은 추가 모형·소광·쌍성 가정이 있을 때 별도 수행합니다.', ['daophot', 'cmd'], True),
    ]))
    result.append(recipe('icl', 'ICL · 은하단 확산광',
        '확산광 보존 이력이 확인된 영상, PSF·마스크·은하단 중심·적색편이 및 대조장.',
        '정의를 명시한 BCG+확산광/ICL 프로파일, 오차·마스크·검증 결과', [
        node('center', '중심·측정 정의', 'BCG 보존/제외, 반경, 단위와 ICL 정의를 고정합니다.', ['icl', 'center']),
        node('mask', '광원·밝은 별 마스크', 'BCG와 위성 처리 규칙을 구별하고 검출 누락을 평가합니다.', ['icl', 'mask']),
        node('sky', '하늘 배경·대조장', '은하단 외곽을 빈 하늘로 가정하지 않습니다. 대조장과 배경 전달을 검증합니다.', ['icl', 'background']),
        node('profile', '표면밝기 프로파일', 'PSF·소광·밴드와 면적을 확인합니다.', ['icl', 'profile']),
        node('stack', '은하단별 프로파일 결합', '같은 은하단 가중, 참여 ID, 반경별 면적을 기록합니다. 스택 입력 manifest를 먼저 준비합니다.'),
        node('validate', 'PSF·주입 회수·공분산', '배경/마스크 변화와 은하단 단위 공분산을 검증합니다.'),
    ]))
    result.append(recipe('lsbg', 'LSBG · 탐색과 구조',
        '넓은 시야 또는 대상 영상, PSF·영점·마스크와 sky 처리 이력.',
        'LSBG 후보 표, 각크기·표면밝기·오차, 컷아웃 및 선택 완전도', [
        node('mask', '밝은 광원·결함 처리', '확장광원을 보존하며 별 날개와 결함을 가립니다.', ['lsbg', 'mask']),
        node('sky', '배경 모형 설정·실행', '모형이 LSBG를 지우는 효과를 인공은하로 확인합니다.', ['lsbg', 'clean']),
        node('detect', 'LSBG 후보 검출', '설정한 배경 처리와 크기 범위를 확인하고 검출합니다.', ['lsbg', 'detect']),
        node('phot', '원본 영상 측광', '검출용 평활 영상과 측광용 원본을 구별합니다.', ['lsbg', 'photometry']),
        node('shape', 'PSF·이웃을 포함한 구조', 'Sérsic 적합과 잔차·퇴화·실패를 검토합니다.', ['lsbg', 'sersic']),
        node('select', '표면밝기·크기·품질 선택', '거리 없이 물리 크기를 확정하지 않습니다.', ['lsbg', 'filter']),
        node('complete', '인공은하·오검출 검증', '검출부터 최종 선택까지 회수하고 순도는 독립 근거로 평가합니다.', ['completeness', 'measure']),
    ]))
    common_spec = [
        node('link', '스펙트럼·대상 연결', '1D/2D/IFU와 파장·flux·오차 단위를 설정하고 대상 catalog와 연결합니다.', ['spectra', 'link']),
        node('calibrate', '추출·파장·분해능·하늘 검토', '2D/IFU는 추출 영역과 sky, 1D는 기존 교정을 확인합니다. 파장 air/vacuum, 관측/휴지계와 중복 보정을 확인합니다.', ['spectra', 'viewer']),
    ]
    result.append(recipe('stellar_spectrum', 'Stellar spectrum · 항성 분광',
        '1D 스펙트럼 또는 교정된 2D/IFU, 파장·flux·오차와 분해능; 대응 대상 목록.',
        '흡수선·시선속도 등 검증된 측정, 스펙트럼·잔차 및 오차', common_spec + [
        node('continuum', '연속광·흡수선 영역', 'telluric·결함 영역, 연속광과 기기 분해능을 확인합니다.'),
        node('science', '선 측정·시선속도', '원하는 측정에 맞는 설정·템플릿과 기준계를 확인합니다.', ['spectra', 'science']),
        node('type', '분광형·대기변수 확장', '온도 letter 도구에는 검증된 사용자 모형이 필요합니다. Teff/logg/금속도 정밀 적합을 대신하지 않습니다.', ['spectra', 'startype'], True),
        node('validate', '잔차·오차·표준별 검증', '중심속도 보정은 시간·관측소 정보가 있을 때만 수행하며 회수 실험·대조 스펙트럼을 사용합니다.'),
    ]))
    result.append(recipe('agn_spectrum', 'AGN spectrum · 활동은하핵 분광',
        '교정된 스펙트럼과 오차·분해능·적색편이 근거, 관측 파장 범위.',
        '연속광·방출선 측정, 분해 모형과 불확실성; 지원되는 범위의 진단', common_spec + [
        node('redshift', '선 후보·적색편이', '넓은 선만으로 정한 적색편이를 좁은 선 등으로 대조합니다.', ['spectra', 'fit']),
        node('decompose', 'AGN·host·Fe II 연속광 분해', '적절한 템플릿 및 넓은/좁은 선 공동 모형이 필요합니다. 전용 adapter는 미구현입니다.'),
        node('lines', '넓은선·좁은선 공동 적합', '블렌딩·기기 분해능·흡수와 공분산을 포함합니다. 단순 Gaussian 결과를 전용 분해의 대체로 표시하지 않습니다.'),
        node('diagnostic', '선비·광도 진단', '필요한 선의 검출/상한과 소광을 전파합니다. 블랙홀 질량은 검증된 보정관계와 광도가 있을 때만 별도 산출합니다.'),
        node('validate', '모형·잔차·오차 검증', '대체 연속광/성분 수와 회수 실험으로 체계오차를 확인합니다.'),
    ]))
    result.append(recipe('galaxy_spectrum', 'Galaxy spectrum · 은하 분광',
        '교정된 은하 스펙트럼 또는 IFU와 오차·분해능·추출 구경 정보.',
        '적색편이·선 측정·오차 및 지원되는 종족/운동학 진단', common_spec + [
        node('redshift', '적색편이·선 후보', '여러 선과 흡수 특징의 일관성을 확인합니다.', ['spectra', 'fit']),
        node('continuum', '별 연속광·흡수 보정', '검증된 항성종족 템플릿이 필요한 분기는 별도 준비합니다.'),
        node('science', '흡수·방출선 측정', 'SFR은 flux 교정·거리·소광·AGN 오염 및 구경 조건을 확인한 경우에만 해석합니다.', ['spectra', 'science']),
        node('kin', '운동학 확장', '2D/IFU와 기기 분해능이 있을 때 속도/분산 지도를 검토합니다.', ['spectra', 'kin'], True),
        node('validate', '잔차·퇴화·오차 검토', '별종족 정밀 추론과 간단한 선 기반 분류를 구별합니다.'),
    ]))
    result.append(recipe('galaxy_morphology', 'Single galaxy · 형태와 구조',
        '대상 은하 영상, PSF·오차·마스크·영점; 물리 크기는 거리 필요.',
        '형태 지표, 광도 프로파일·모형·잔차, 오차 및 해석 범위', [
        node('detect', '대상·이웃 검출', '대상 은하를 선택하고 주변 천체를 식별합니다.', ['extract', 'extract']),
        node('psf', 'PSF·배경·마스크', '대상 외곽을 보호하며 sky/PSF 모형과 이웃 처리를 확인합니다.', ['star_psf', 'build_psf']),
        node('shape', 'Sérsic 구조 측정', '반경·축비·위치각과 잔차를 검토하고 PSF 처리 적용 여부를 확인합니다.', ['morphology', 'sersic']),
        node('nonparam', 'CAS/Gini/M20', '분할 영역·중심·배경과 신호대잡음 의존성을 점검합니다.', ['morphology', 'morphometry']),
        node('validate', '모형·배경 민감도', '해상도와 주입 실험으로 측정 편향을 평가합니다. 다성분 분해는 지원·검증 범위에 따라 별도 추가합니다.'),
    ]))
    result.append(recipe('source_extraction', 'Source extractor · 광원 검출',
        '영상, mask/variance, pixel scale·영점과 검출 목적.',
        '검출 목록·segmentation, 측광/위치 오차·플래그·완전도', [
        node('configure', '배경·threshold·deblend 설정', '점광원/확장광원 목적과 잡음·필터를 맞춥니다.', ['settings', 'extract']),
        node('detect', '광원 검출·측광', '설정된 backend를 기록합니다. SEP 기반 추출과 외부 SExtractor 실행을 구별합니다.', ['extract', 'extract']),
        node('segment', '분할 영역·블렌딩 검토', '분할/병합·밝은 별·영상 결함을 확인합니다.', ['photometry', 'segmap']),
        node('complete', '검출 성능·한계', '목적에 맞는 점광원/확장광원을 넣어 검출 및 품질 컷의 회수율을 측정합니다.', ['completeness', 'measure']),
    ]))
    epochs = [
        node('epochs', '시점 영상·시간 설정', '노출별 시간계·관측소·필터·WCS·오차를 기록합니다.', ['moving', 'select']),
        node('align', '측성·광도 정렬', '정렬 오차·seeing·영점과 유효 영역을 검토합니다.', ['moving', 'align']),
        node('difference', '기준영상·차영상', '기준영상이 새 천체 신호를 포함하는지 확인합니다. PSF와 분산을 전파합니다.', ['moving', 'difference']),
        node('detect', '시점별 후보 검출', '우주선·결함·회절선·쌍극 잔차를 제거합니다.', ['moving', 'detect']),
    ]
    result.append(recipe('asteroid', '소행성 · 이동천체 탐색',
        '여러 시점의 영상과 정확한 관측시각·관측소; 기존 영상 목록을 불러옵니다.',
        '이동 후보·tracklet·측성/측광 표 및 가능한 범위의 궤도 제약', epochs + [
        node('link', '운동 일관성·tracklet', '시간순 위치와 속도 범위·오연결률을 검토합니다.', ['moving', 'link']),
        node('identify', '알려진 천체 대응', '관측시각·위치 오차를 포함해 알려진 천체를 대조합니다.', ['moving', 'identify']),
        node('orbit', '궤도·불확실성', '짧은 관측호는 유일 궤도를 보장하지 않습니다.', ['moving', 'orbit'], True),
        node('validate', '후보 검토·주입 회수', '이동 인공천체의 회수와 오연결을 확인합니다. 외부 관측 보고는 별도 사용자 작업입니다.'),
    ]))
    result.append(recipe('transient', 'Transient · 일시천체 탐색',
        '새 영상과 적절한 기준영상, 다중 시점 자료·시간·PSF·오차.',
        '후보 목록, 차영상·광도곡선·비검출 상한 및 검증 기록', epochs + [
        node('candidates', '고정 위치 후보·이동체 배제', '양/음 잔차, 위치 지속성·이동체·알려진 변광을 검토합니다.', ['moving', 'transients']),
        node('lightcurve', '강제 측광·광도곡선', '비검출과 상한, 기준영상 flux·공통 오차를 포함합니다.', ['moving', 'lightcurve']),
        node('classify', '후보 분류 검토', '분류기는 훈련 분포·관측 조건을 검증해야 합니다. 후보 등급을 확정 분류로 바꾸지 않습니다.', ['lightcurves', 'classify'], True),
        node('validate', '반복 관측·오검출·회수', '인공 transient 회수와 독립 검토를 수행합니다. 알림·보고 전송은 별도 사용자 작업입니다.'),
    ]))
    agn = next(item for item in result if item['id'] == 'agn_spectrum')
    inference = node('inference', '분류·소광·광도·블랙홀 질량 통합 해석', '선행 공동 적합을 직접 연결합니다. inference-config JSON에 좁은선 BPT 선택, 소광법칙 k/외부 E(B-V) 또는 Balmer 비, 광도거리·물리 flux 변환, 확인된 넓은 Hα/기기 FWHM/Greene-Ho 2005 식6 적용근거·scatter를 입력합니다. 전체 적합 공분산과 독립 보조입력 오차를 전파하고 각 항목의 결과/모호성/계산불가 사유를 저장합니다. 빈 설정은 미계산이며 결과는 조건부 추론입니다.', ['ds10core', 'spectrum-inference'])
    inference.update(web_action=['ds10_spectrum', 'inference'], web_detail=inference['detail'], input_of='lines', input_artifact='joint.json', desktop_input_param='spectrum-inference-input')
    agn['steps'].insert(next(i for i,step in enumerate(agn['steps']) if step['id'] == 'validate'), inference)
    native = {
        ('extract', 'extract'): ['ds10_detect', 'detect'],
        ('morphology', 'morphometry'): ['ds10_morph', 'morph'],
        ('star_psf', 'build_psf'): ['ds10_psfstamp', 'psfstamp'],
        ('daophot', 'phot'): ['ds10_aper_phot', 'aper-phot'],
    }
    for item in result:
        if item['id'].endswith('_spectrum'):
            target = 'science' if item['id'] == 'stellar_spectrum' else 'redshift'
            position = next(i for i, step in enumerate(item['steps']) if step['id'] == target)
            proposal = node('search', '선 탐색·적색편이 후보 목록', '초기 z 없이 기준선 목록과 관측 피크를 대조합니다. 최소 두 선의 후보·대응표와 대안을 저장합니다. 점수는 확률이 아니며, 후보 검토 후 suggested_lines와 z를 다음 정밀 선 측정 설정에 입력하세요. 자동 확정·기존 설정 변경은 하지 않습니다.', ['ds10core', 'spectrum-search'])
            proposal['web_action'] = ['ds10_spectrum', 'search']
            proposal['web_detail'] = proposal['detail']
            item['steps'].insert(position, proposal)
        for step in item['steps']:
            if item['id'].endswith('_spectrum'):
                if step['id'] == 'load':
                    step['web_action'] = ['spectrum_load']
                    step['web_import_action'] = ['ds10_spectrum', 'import']
                    step['detail'] += ' 교정된 단일 1D TSV: WAVE_A/FLUX/FLUX_ERR 및 선택 MASK(0=사용, 1=제외)를 지원합니다. 2D/IFU는 먼저 별도 추출·교정합니다.'
                if step['id'] == 'audit':
                    step['action'] = ['ds10core', 'spectrum-audit']
                    step['web_action'] = ['ds10_spectrum', 'audit']
                    step['detail'] = '1D TSV의 증가하는 Angstrom 파장, 유효 flux·양수 오차·마스크를 검사합니다. air/vacuum과 flux 단위를 명시해야 하며 교정·추출은 수행하지 않습니다.'
                if step['id'] == 'continuum':
                    step['title'] = '1D 연속광·잔차·정규화'
                    step['action'] = ['ds10core', 'spectrum-continuum']
                    step['web_action'] = ['ds10_spectrum', 'continuum']
                    step['detail'] = '제외 파장 구간·오차 가중·반복 clipping으로 Chebyshev 연속광을 적합합니다. 연속광·조건부 모형 공분산·잔차·정규화 스펙트럼과 원본을 저장합니다. 항성종족/host/Fe II 분해 및 흡수 보정을 대신하지 않습니다.'
                if step['id'] in ('science', 'redshift'):
                    task = 'redshift' if step['id'] == 'redshift' else 'lines'
                    step['action'] = ['ds10core', 'spectrum-' + task]
                    step['web_action'] = ['ds10_spectrum', task]
                    step['title'] = '동정된 여러 선의 적색편이 후보' if task == 'redshift' else '동정된 고립선 측정'
                    step['detail'] = '사용자가 NAME:REST_A로 동정한 고립선과 초기 적색편이를 입력합니다. 관측 반창 폭(Angstrom), 방출/흡수와 최소 S/N을 설정합니다. 국소 선형 연속광+Gaussian을 공동 적합하고 flux·중심·관측 폭·공분산·잔차를 저장합니다. 두 개 이상 승인된 선이 일치할 때만 조건부 적색편이 후보를 제공합니다. 자동 선 동정/전역 적색편이 검색·시선속도 보정·AGN 분해는 수행하지 않습니다. flux는 Angstrom당 밀도 단위여야 합니다.'
                    step['web_detail'] = step['detail']
                if step['id'] in ('type', 'kin'):
                    step['web_action'] = []
                    step['web_detail'] = '교정된 1D TSV 경로의 전용 과학 측정 adapter는 아직 미구현입니다. Desktop 대상 catalog/FITS용 도구에 이 표를 바로 전달하지 않습니다. 현재 입력 검사·연속광 결과를 저장해 후속 분석에 사용하세요.'
                if step['id'] in ('link', 'calibrate'):
                    step['web_action'] = ['settings', 'ds10_spectrum']
                    step['web_result_of'] = 'audit'
                    step['web_detail'] = '웹 1D TSV 경로의 단위·파장계·입력 검사 산출물 검토입니다. 영상 대상 연결 및 2D/IFU 추출·sky 교정은 별도 준비해야 합니다.'
            if item['id'] == 'agn_spectrum' and step['id'] == 'decompose':
                step['title'] = '사용자 host·Fe II 템플릿 세기 분해'
                step['action'] = ['ds10core', 'spectrum-decompose']
                step['web_action'] = ['ds10_spectrum', 'decompose']
                step['detail'] = '입력 TSV에 같은 관측 파장 격자의 TEMPLATE_HOST/TEMPLATE_FEII를 포함하고 출처·파장계·적색편이·분해능 준비를 template-note에 기록합니다. 고정 F_lambda power-law index/pivot와 사용자 템플릿의 양수 세기를 적합합니다. 방출선을 exclude로 제외하고 fit-range를 지정하세요. 모형·조건부 공분산·잔차·원본을 저장하며 템플릿 모양/폭/소광은 적합하지 않습니다. 다음 공동 선 적합의 continuum-model=templates에서 세기를 다시 함께 적합해 선과의 공분산을 반영하세요.'
                step['web_detail'] = step['detail']
            if item['id'] == 'agn_spectrum' and step['id'] == 'diagnostic':
                step['title'] = '선비·log 비와 공분산'
                step['action'] = ['ds10core', 'spectrum-diagnostic']
                step['web_action'] = ['ds10_spectrum', 'diagnostic']
                step['input_of'] = 'lines'
                step['input_artifact'] = 'joint.json'
                step['desktop_input_param'] = 'spectrum-diagnostic-input'
                step['detail'] = '앞선 공동 적합 JSON을 자동 연결합니다. ratios에 LABEL=성분+성분/성분+성분을 입력하면 합계 선비·log10 비와 전체 성분 공분산을 전파한 오차/선비 간 공분산을 저장합니다. flagged 성분·낮은 합계 S/N은 계산하지 않습니다. 관측 선비이며 소광/광도/분류/블랙홀 질량은 산출하지 않습니다.'
                step['web_detail'] = step['detail']
            if item['id'] == 'agn_spectrum' and step['id'] == 'lines':
                step['title'] = '방출선 성분·국소 연속광 공동 적합'
                step['action'] = ['ds10core', 'spectrum-joint']
                step['web_action'] = ['ds10_spectrum', 'joint']
                step['detail'] = '사용자 NAME:REST_A:SIGMA_MIN_A:SIGMA_MAX_A 성분과 관측 fit-range를 입력합니다. 선 중심은 공통 적색편이에 묶고 국소 선형 연속광+여러 방출 Gaussian을 공동 적합합니다. 좁은/넓은 성분은 폭 범위를 지정하며 같은 기준 파장에서는 서로 겹치지 않는 폭 범위를 사용합니다. 전체/성분 flux 공분산·잔차·마스크·원본을 저장합니다. continuum-model=templates이면 같은 입력의 사용자 host/Fe II 템플릿+고정 power-law 세기를 선과 함께 재적합합니다. 템플릿 모양/분해능은 사전 준비하며 기기 폭 제거·블랙홀 질량은 지원하지 않습니다.'
                step['web_detail'] = step['detail']
            alternate = native.get(tuple(step['action']))
            if alternate:
                step['web_action'] = alternate
            if item['id'] == 'source_extraction' and step['id'] == 'configure':
                step['web_action'] = ['settings', 'ds10_detect']
            if item['id'] == 'source_extraction' and step['id'] == 'segment':
                step['web_action'] = ['settings', 'ds10_detect']
                step['web_result_of'] = 'detect'
                step['web_detail'] = '검출 작업의 segmap.fits와 catalog.tsv를 검토합니다. 분할 영역은 동일 검출에서 생성되며 별도 재검출하지 않습니다.'
            if item['id'] == 'source_extraction' and step['id'] == 'complete':
                step['web_action'] = ['ds10_counts', 'counts']
                step['web_detail'] = '내장 검출기의 인공 점광원/은하 주입 회수입니다. 검출 단계의 threshold·평활·배경·deblend와 이 도구의 설정을 일치시키세요. PSF 측광이나 통합 LSBG 선택의 완전도를 대신하지 않습니다.'

            if item['id'] in ('globular_cluster', 'open_cluster'):
                if step['id'] in ('fit', 'phot'):
                    step['web_action'] = ['ds10_psf_phot', 'psf-phot']
                    step['web_detail'] = '웹 통합 PSF 측광: 별 선택·ePSF 구성·이웃 공동 적합·구경 보정을 한 번에 실행합니다. 잔차, 보정 플래그와 영점 근거를 확인하세요.'
                if step['id'] in ('pickpsf', 'psf', 'calibrate'):
                    step['web_action'] = ['settings', 'ds10_psf_phot']
                    step['web_detail'] = '웹에서는 통합 PSF 측광 단계에서 계산합니다. 이 단계는 해당 설정과 산출물 검토이며 독립 계산이 아닙니다.'
                    step['web_result_of'] = 'fit' if item['id'] == 'globular_cluster' else 'phot'
            if item['id'] == 'lsbg':
                if step['id'] == 'detect':
                    step['web_action'] = ['ds10_lsbg', 'find']
                    step['web_detail'] = '웹 통합 LSBG 분석: 밝은 광원 마스크·배경·검출·측광·구조 추정·후보 선택을 한 번 실행합니다. 주입 회수나 PSF 포함 정밀 모형의 검증을 대신하지 않습니다.'
                if step['id'] in ('mask', 'sky', 'phot', 'shape', 'select'):
                    step['web_action'] = ['settings', 'ds10_lsbg']
                    step['web_result_of'] = 'detect'
                    step['web_detail'] = '웹에서는 통합 LSBG 분석의 설정·산출물을 검토합니다. 이 단계의 별도 실행은 없으며 구조 추정의 PSF 처리와 편향을 따로 확인해야 합니다.'
            if item['id'] == 'icl' and step['id'] == 'stack':
                step['action'] = ['ds10core', 'profile-stack']
                step['web_action'] = ['ds10_workflow', 'stack']
                step['detail'] += ' 공통 kpc 격자의 CLUSTER_ID/R_KPC/SB/SB_ERR/AREA_KPC2 및 SB_UNIT/BAND/ICL_DEFINITION/SKY_METHOD/PSF_TARGET TSV가 필요합니다. 은하단 단위 bootstrap 공분산과 참여 ID를 저장합니다.'
    return result


def validate_catalog(items, root=None):
    ids = set()
    for r in items:
        if r['id'] in ids:
            raise ValueError('duplicate recipe')
        ids.add(r['id'])
        nodes = [s['id'] for s in r['steps']]
        if len(nodes) != len(set(nodes)):
            raise ValueError('duplicate step')
        for a, b in r['edges']:
            if a not in nodes or b not in nodes or nodes.index(a) >= nodes.index(b):
                raise ValueError('invalid edge')
        if root:
            for s in r['steps']:
                a = s['action']
                if len(a) == 2:
                    pid = a[1] if a[0] == 'settings' else a[0]
                    p = Path(root) / 'plugins' / pid / 'plugin.json'
                    m = json.loads(p.read_text())
                    if a[0] != 'settings' and a[1] not in {x['id'] for x in m['steps']}:
                        raise ValueError(f'unknown action {a}')


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def checkpoint(destination, state, outputs=(), inputs=()):
    """Publish a new snapshot only after all explicit outputs were copied and hashed."""
    destination = Path(destination).absolute()
    destination.mkdir(parents=True, exist_ok=True)
    # Exclusive writer; an interrupted writer leaves evidence, never a valid snapshot.
    lock = destination / '.saving'
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    os.close(fd)
    staging = None
    try:
        files = list(dict.fromkeys(str(Path(p).resolve(strict=True)) for p in outputs))
        for name in files:
            if not Path(name).is_file():
                raise ValueError(f'output must be an explicit regular file: {name}')
        total = sum(Path(p).stat().st_size for p in files)
        if socket.gethostname().split('.')[0] == 'syntax' and total > 256 * 1024 * 1024:
            raise ValueError('Large checkpoint copying must run on lageunha, not the syntax login node')
        staging = Path(tempfile.mkdtemp(prefix='.pending-', dir=destination))
        (staging / 'artifacts').mkdir()
        inventory = []
        for i, name in enumerate(files):
            src = Path(name)
            before = src.stat()
            relative = f'artifacts/{i:04d}_{src.name}'
            target = staging / relative
            shutil.copyfile(src, target)
            sha = digest(target)
            after = src.stat()
            if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino) or digest(src) != sha:
                raise ValueError(f'output changed while saving: {name}')
            inventory.append(dict(source=name, path=relative, bytes=after.st_size, sha256=sha))
        references = []
        for name in dict.fromkeys(inputs):
            p = Path(name).resolve(strict=True)
            if not p.is_file():
                raise ValueError(f'input must be a file: {p}')
            st = p.stat()
            references.append(dict(path=str(p), bytes=st.st_size, mtime_ns=st.st_mtime_ns,
                                   policy='external_reference_size_mtime_only'))
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
        snapshot = dict(schema_version=1, created=stamp, state=state,
                        artifacts=inventory, inputs=references,
                        science_validated=False, snapshot_complete=True)
        manifest = staging / 'session.json'
        with manifest.open('w', encoding='utf-8') as f:
            json.dump(snapshot, f, ensure_ascii=False, indent=2)
            f.write('\n')
            f.flush()
            os.fsync(f.fileno())
        final = destination / f'snapshot-{stamp}-{uuid.uuid4().hex[:8]}'
        staging.rename(final)
        staging = None
        return final / 'session.json'
    finally:
        if staging:
            shutil.rmtree(staging)
        lock.unlink()


def restore(manifest, work=None):
    manifest = Path(manifest).resolve(strict=True)
    data = json.loads(manifest.read_text(encoding='utf-8'))
    if data.get('schema_version') != 1 or data.get('snapshot_complete') is not True:
        raise ValueError('unsupported or incomplete snapshot')
    if socket.gethostname().split('.')[0] == 'syntax' and sum(a['bytes'] for a in data['artifacts']) > 256 * 1024 * 1024:
        raise ValueError('Large checkpoint verification must run on lageunha, not the syntax login node')
    for item in data['artifacts']:
        p = (manifest.parent / item['path']).resolve(strict=True)
        if not p.is_relative_to(manifest.parent) or not p.is_file():
            raise ValueError('invalid artifact path')
        if p.stat().st_size != item['bytes'] or digest(p) != item['sha256']:
            raise ValueError(f'checkpoint artifact damaged: {item["path"]}')
        item['resolved_path'] = str(p)
    changes = []
    for ref in data['inputs']:
        p = Path(ref['path'])
        if not p.is_file():
            changes.append(dict(path=str(p), reason='missing'))
        else:
            st = p.stat()
            if st.st_size != ref['bytes'] or st.st_mtime_ns != ref['mtime_ns']:
                changes.append(dict(path=str(p), reason='changed'))
    data['input_changes'] = changes
    if work is not None:
        work = Path(work).absolute()
        work.mkdir(parents=True, exist_ok=False)
        try:
            for i, item in enumerate(data['artifacts']):
                target = work / f'{i:04d}_{Path(item["path"]).name}'
                shutil.copyfile(item['resolved_path'], target)
                if digest(target) != item['sha256']:
                    raise ValueError('artifact changed during restore')
                item['resolved_path'] = str(target)
        except Exception:
            shutil.rmtree(work)
            raise
    return data
