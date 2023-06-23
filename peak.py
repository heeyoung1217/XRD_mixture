import numpy as np

# 주어진 데이터
x_data = np.array([18.31, 30.12, 35.47, 37.11, 43.11, 47.2, 53.49, 57.02, 62.61, 65.83, 71.03, 74.07,
                   75.08, 79.04, 81.98, 86.84, 89.75])
y_data = np.array([8.99977091, 26.87434724, 100., 4.4026637, 27.14741223, 0.82698172,
                   9.57221282, 33.12549803, 50.76984533, 1.2826262, 3.74448687,
                   9.89741833, 3.12190021, 4.05311885, 0.69598005, 4.34431795, 16.99921404])

# x 데이터와 y 데이터 매칭
data = list(zip(x_data, y_data))

# 새로운 x 데이터 생성
new_x_data = np.arange(10.00, 80.01, 0.01)

# y 데이터 대입
new_y_data = []
for x in new_x_data:
    closest_data = min(data, key=lambda d: abs(d[0] - x))
    new_y_data.append(closest_data[1])

# 결과 출력
for x, y in zip(new_x_data, new_y_data):
    print("{:.2f}: {:.8f}".format(x, y))