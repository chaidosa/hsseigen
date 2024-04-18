
para = """
    Completed Hss tree
Reached superDC
Success Divide

Done. 0.002336 usecs
Completed Hss tree
Reached superDC
Success Divide

Done. 0.000928 usecs
Completed Hss tree
Reached superDC
Success Divide

Done. 0.001298 usecs
Completed Hss tree
Reached superDC
Success Divide

Done. 0.001304 usecs
Completed Hss tree
Reached superDC
Success Divide

Done. 0.001944 usecs
"""

import re

def extract_floats(text):
    # Regular expression to find float numbers
    float_pattern = r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
    floats = re.findall(float_pattern, text)
    return [float(num) for num in floats]

def calculate_average(numbers):
    if not numbers:
        return None
    return sum(numbers) / len(numbers)

def main():
    paragraph = """
    Lorem ipsum dolor sit amet, consectetur adipiscing elit. 
    Sed at 2.5 nulla nec sapien volutpat condimentum. 
    Cras 1.5 eget justo at libero 3.7 tempor convallis. 
    Phasellus 4.2 vel quam nec justo pretium ultricies.
    """

    floats = extract_floats(para)
    average = calculate_average(floats)

    if average is not None:
        print("Extracted float numbers:", floats)
        print("Average:", average)
        print("Sum:", sum(floats))
    else:
        print("No float numbers found in the text.")

if __name__ == "__main__":
    main()
