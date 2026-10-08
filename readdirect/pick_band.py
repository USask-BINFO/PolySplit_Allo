import sys

CI = int(sys.argv[2]) if len(sys.argv) > 2 else 3
CAP = int(sys.argv[3]) if len(sys.argv) > 3 else 300
d = {}
for ln in open(sys.argv[1]):
    p = ln.split()
    if len(p) >= 2:
        try:
            c, n = int(p[0]), int(p[1])
        except ValueError:
            continue
        if CI <= c <= CAP:
            d[c] = n
cs = sorted(d)
if len(cs) < 4:
    print(CI, 3 * CI); sys.exit()

# valley = first local minimum after the error-tail decline
valley = None
for i in range(1, len(cs) - 1):
    if d[cs[i - 1]] > d[cs[i]] and d[cs[i + 1]] >= d[cs[i]]:
        valley = cs[i]; break

if valley is None:
    flatten = cs[0]
    for i in range(1, len(cs)):
        if d[cs[i]] > 0.97 * d[cs[i - 1]]:
            flatten = cs[i]; break
    print(int(flatten), int(min(3 * flatten, CAP))); sys.exit()

window = [c for c in cs if valley < c <= valley + max(10, valley)]
peak = max(window, key=lambda c: d[c]) if window else valley + max(3, valley // 4)
print(int(valley), int(min(2 * peak - valley, CAP)))
